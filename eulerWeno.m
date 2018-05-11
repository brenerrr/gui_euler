function [ dF_dx ] = eulerWeno( inp, U, F )
% Calculates the derivative of 1-D euler fluxes with a classic weno scheme

if inp.WENO==1
    f_p=zeros(3,inp.SIZEX);
    f_n=f_p;
    T_inv = zeros(3,3,inp.SIZEX);
    T=T_inv;
    G=zeros(3,inp.SIZEX);
    W=G;
    W_p=G;
    W_n=G;
    w_p=G;
    w_n=G;
    
    u = U(2,:)./U(1,:);
    P = (inp.GAMA - 1) * (U(3,:) - 0.5*U(1,:).*u.^2);    
    c= sqrt(inp.GAMA*P/U(1,:));    
    lambda(:,:) = [ u; u+c; u-c];
    coef = max(abs(lambda'));

    for i=inp.V1-1:inp.VN    
         
        U_temp = (U(:,i)+U(:,i+1))/2; % properties must be evaluated at borders
        
        aux_temp(1) = U_temp(2) / U_temp(1); % velocity
        aux_temp(2) = (inp.GAMA - 1) * (U_temp(3) - 0.5*U_temp(1).*aux_temp(1).^2); % pressure
        c_temp = sqrt(inp.GAMA*aux_temp(2)/U_temp(1)); % sound's speed
        
        [R, R_inv]=jacob( U_temp(1), aux_temp(1), c_temp, inp.GAMA);

        % Decouple U and F 
        T(:,:,i)=R(:,:);
        G(:,i-2:i+3)=R_inv * U(:,i-2:i+3);
        W(:,i-2:i+3)=R_inv * F(:,i-2:i+3); 
        
        % Get positive and negative fluxes
        for j=i-2:i+3
            W_p(:,j) = 0.5*(W(:,j) + coef(:).*G(:,j) );
            W_n(:,j) = 0.5*(W(:,j) - coef(:).*G(:,j) );         
        end
        
        % Apply weno scheme in each flux
        for j=1:3
            q = inp.N_N_POINTS;
            w_p(j,i) = fluxWenoPos(W_p(j,i-q:i+q),inp);
            w_n(j,i) = fluxWenoNeg(W_n(j,i-q+1:i+q+1),inp);
        end
      
    end

    % Couple everything back up
    for i=inp.V1-1:inp.VN
       R(:,:) = T(:,:,i);
       f_p(:,i) = R * w_p(:,i);   %fluxo positivo
       f_n(:,i) = R * w_n(:,i);   %fluxo negativo
    end
    
	% Calculate final flux
    flux = f_p(:,inp.V1-1:inp.VN) + f_n(:,inp.V1-1:inp.VN) ;

end

dF_dx=( flux(:,2:end) - flux(:,1:end-1) ) ./ inp.DELTA_X;

end

%% Subroutines

function [R, R_inv, e_val] = jacob(rho, u, c, gama)
% Calculates the transformation matrix R, R_inv and the eigenvalues of a
% 1-D euler jacobian

alpha=(u^2) * 0.5;
beta=gama-1;

S = [1 0 0 ; ...
    -u/rho 1/rho 0 ; ...
    alpha*beta -u*beta beta];

Ca = [1 0 -1/c^2 ; ...
      0 rho*c 1 ; ...
      0 -rho*c 1];

S_inv = [1 0 0 ; ...
         u rho 0 ; ...
         alpha rho*u 1/beta];
     
Ca_inv = [1 0.5*1/c^2 0.5*1/c^2 ; ...
          0 0.5*1/(rho*c) -0.5*1/(rho*c) ; ...
          0 0.5 0.5];

                   
e_val=zeros(1,3);
e_val(1) = u;
e_val(2) = u + c;
e_val(3) = u - c;

R = S_inv * Ca_inv;
R_inv = Ca * S;

end

function [vhr] =fluxWenoPos (p, inp)
% Calculate positive fluxes of p at i+1/2
vhr=0;
vhr_dummy=zeros(inp.TOTAL_POINTS,1);

    if isreal(p)==0
        error ('imaginario') 
    end
    
    for l=1:inp.TOTAL_POINTS
        
        r = inp.TOTAL_POINTS - l + 2;
        vhr_dummy(l) = inp.CRJ(r,:)*p(l:l+inp.N_N_POINTS)';
                    
    end
    
    beta=smoothFactor(p, inp.TOTAL_POINTS);
    
    if inp.TOTAL_POINTS == 2
       C = [1/3 ; 2/3]; 
    end
    
    if inp.TOTAL_POINTS == 3
       C = [1/10 ; 3/5 ; 3/10]; 
    end
    alpha=C./((inp.EPSILON+beta).^2);
    
    w=alpha./sum(alpha);
    
    
    for m=1:inp.TOTAL_POINTS
        vhr=vhr+w(m)*vhr_dummy(m);
        if w(m)<0
           error ('peso negativo') 
        end
    end
    
end

function [vhl] =fluxWenoNeg (p, inp)
% Calculate negative fluxes of p at i+1/2

vhl=0;

vhl_dummy=zeros(inp.TOTAL_POINTS,1);

    if isreal(p)==0
        error ('imaginario') 
    end
    
    for l=1:inp.TOTAL_POINTS
        r = inp.TOTAL_POINTS - l + 2;
        vhl_dummy(l) = inp.CRJ(r-1,:)*p(l:l+inp.N_N_POINTS)';
    end
    
    beta=smoothFactor(p, inp.TOTAL_POINTS);
    
    if inp.TOTAL_POINTS == 2
       C_l = [2/3 ; 1/3];
    end
    
    if inp.TOTAL_POINTS == 3
       C_l = [3/10 ; 3/5 ; 1/10]; 
    end
    alpha_l=C_l./((inp.EPSILON+beta).^2);
    
    w_l=alpha_l./sum(alpha_l);
    
    
    for m=1:inp.TOTAL_POINTS
        vhl=vhl+w_l(m)*vhl_dummy(m);
        if w_l(m)<0
           error ('peso negativo') 
        end
    end
    

end

function [beta]=smoothFactor(f, s)
% Calculates smooth factor for second or forth order weno interpolation

    if s == 2
        beta=[ (f(2)-f(1))^2 ; (f(3)-f(2))^2 ];
    end
        
    if s == 3
        beta = [ 13/12*(f(1)-2*f(2)+f(3))^2+1/4*(f(1)-4*f(2)+3*f(3))^2 ;...
                 13/12*(f(2)-2*f(3)+f(4))^2+1/4*(f(2)-f(4))^2 ;...
                 13/12*(f(3)-2*f(4)+f(5))^2+1/4*(3*f(3)-4*f(4)+f(5))^2];
    end

end
