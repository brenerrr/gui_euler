function [ dF_dx ] = Roe( inp, U, F )
% Calculates the derivative of the 1-D euler fluxes using a roe method that
% can be put together with weno interpolation

aux(1,1:inp.SIZEX) = U(2,:) ./ U(1,:); % velocidade
aux(2,1:inp.SIZEX) = (inp.GAMA - 1) * (U(3,:) - 0.5*U(1,:).*aux(1,:).^2); % pressao
c = sqrt(inp.GAMA*aux(2,:)./U(1,:));

%%  Roe with WENO
if inp.WENO==1
    
    U_r=zeros(3,inp.SIZEX);
    U_l=U_r;
    aux_r=zeros(1,inp.SIZEX);
    aux_l=aux_r;
    
    [U_r(1,:), U_l(1,:)] = boundariesWeno (U(1,:), inp); 
    [U_r(2,:), U_l(2,:)] = boundariesWeno (U(2,:), inp);
    [U_r(3,:), U_l(3,:)] = boundariesWeno (U(3,:), inp);
    
    aux_r(1,1:inp.SIZEX) = U_r(2,:) ./ U_r(1,:); % velocity
    aux_r(2,1:inp.SIZEX) = (inp.GAMA - 1) * (U_r(3,:) - 0.5*U_r(1,:).*aux_r(1,:).^2); % pressure
    
    aux_l(1,1:inp.SIZEX) = U_l(2,:) ./ U_l(1,:); 
    aux_l(2,1:inp.SIZEX) = (inp.GAMA - 1) * (U_l(3,:) - 0.5*U_l(1,:).*aux_l(1,:).^2);
    
    % Constant boundary conditions
    for i=1:inp.V1-1
       U_r(:,i)=U_r(:,inp.V1);
       U_l(:,i)=U_l(:,inp.V1);
       aux_r(:,i)=aux_r(:,inp.V1);
       aux_l(:,i)=aux_l(:,inp.V1);
       
       U_r(:,inp.VN+i)=U_r(:,inp.VN);
       U_l(:,inp.VN+i)=U_l(:,inp.VN);
       aux_r(:,inp.VN+i)=aux_r(:,inp.VN);
       aux_l(:,inp.VN+i)=aux_l(:,inp.VN);
       
    end
    
    F_r(1,:) = U_r(2,:);
    F_r(2,:) = U_r(2,:) .* aux_r(1,:) + aux_r(2,:);
    F_r(3,:) = ( U_r(3,:) + aux_r(2,:) ) .* aux_r(1, :);
    
    F_l(1,:) = U_l(2,:);
    F_l(2,:) = U_l(2,:) .* aux_l(1,:) + aux_l(2,:);
    F_l(3,:) = ( U_l(3,:) + aux_l(2,:) ) .* aux_l(1, :);
    
    
    % +1/2 Part
    rho_r (inp.V1-1:inp.VN) = ( U_r(1,inp.V1-1:inp.VN) .* U_l(1,inp.V1:inp.VN+1) ).^0.5;
    
    u_r (inp.V1-1:inp.VN) = ( U_r(1,inp.V1-1:inp.VN).^0.5 .* aux_r(1,inp.V1-1:inp.VN) +...
        U_l(1,inp.V1:inp.VN+1).^0.5 .* aux_l(1,inp.V1:inp.VN+1) ) ./ ...
        ( U_r(1,inp.V1-1:inp.VN).^0.5 + U_l(1,inp.V1:inp.VN+1).^0.5 );
    
    h_r (inp.V1-1:inp.VN) = ( U_r(1,inp.V1-1:inp.VN).^0.5 .* ( U_r(3,inp.V1-1:inp.VN) + aux_r(2,inp.V1-1:inp.VN) ) ./ U_r(1,inp.V1-1:inp.VN) + ...
        U_l(1,inp.V1:inp.VN+1).^0.5 .* ( U_l(3,inp.V1:inp.VN+1) + aux_l(2,inp.V1:inp.VN+1) ) ./ U_l(1,inp.V1:inp.VN+1) ) ./ ...
        ( U_r(1,inp.V1-1:inp.VN).^0.5 + U_l(1,inp.V1:inp.VN+1).^0.5 );

else
%% Classic Roe

    % +1/2 Part
    rho_r (inp.V1-1:inp.VN) = ( U(1,inp.V1-1:inp.VN) .* U(1,inp.V1-1+1:inp.VN+1) ).^0.5;
    
    u_r (inp.V1-1:inp.VN) = ( U(1,inp.V1-1:inp.VN).^0.5 .* aux(1,inp.V1-1:inp.VN) +...
        U(1,inp.V1-1+1:inp.VN+1).^0.5 .* aux(1,inp.V1-1+1:inp.VN+1) ) ./ ...
        ( U(1,inp.V1-1:inp.VN).^0.5 + U(1,inp.V1-1+1:inp.VN+1).^0.5 );
    
    h_r (inp.V1-1:inp.VN) = ( U(1,inp.V1-1:inp.VN).^0.5 .* ( U(3,inp.V1-1:inp.VN) + aux(2,inp.V1-1:inp.VN) ) ./ U(1,inp.V1-1:inp.VN) + ...
        U(1,inp.V1-1+1:inp.VN+1).^0.5 .* ( U(3,inp.V1-1+1:inp.VN+1) + aux(2,inp.V1-1+1:inp.VN+1) ) ./ U(1,inp.V1-1+1:inp.VN+1) ) ./ ...
        ( U(1,inp.V1-1:inp.VN).^0.5 + U(1,inp.V1-1+1:inp.VN+1) );
    
end

%% Calculation

    F_rfinal=zeros (3, inp.SIZEX);
    
for i=inp.V1-1:inp.VN
  
    C_r = ( (inp.GAMA-1) * (h_r(i)-0.5*u_r(i)^2) )^0.5;
    
    lambda_r = [aux(1,i+1), aux(1,i+1)+c(i+1), aux(1,i+1)-c(i+1)];
    lambda_l = [aux(1,i), aux(1,i)+c(i), aux(1,i)-c(i)];
    lambda_hat = [u_r(i), u_r(i)+C_r, u_r(i)-C_r];
    
    epsilon=zeros(1,3);
    for j=1:3
        epsilon(j) = 20*max([0, (lambda_hat(j)-lambda_l(j)), (lambda_r(j)-lambda_hat(j)) ]);
    end
    
    A_r=jacob(rho_r(i), u_r(i), C_r, inp.GAMA, epsilon);    
    
    if inp.WENO == 1
        F_rfinal(:,i) = 0.5 * ( F_r(:,i) + F_l(:,i+1) ) - 0.5 * ( A_r * ( U_l(:,i+1) - U_r(:,i) ) );
    else     
        
    F_rfinal(:,i) = 0.5 * ( F(:,i) + F(:,i+1) ) - 0.5 * ( A_r * ( U(:,i+1) - U(:,i) ) );    
    
    end
   
end

dF_dx=( F_rfinal(:,inp.V1:inp.VN) - F_rfinal(:,inp.V1-1:inp.VN-1) ) ./ inp.DELTA_X;
end

%% Functions

function [A] = jacob(rho, u, c, gama, epsilon)
% Calculates jacobian A for an 1-D euler equations

alpha=u^2 * 0.5;
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

for i = 1:3
   if abs(e_val(i)) < epsilon(i)
      e_val(i) =  (e_val(i)^2 + epsilon(i)^2) / (2*epsilon(i)) ;
   end
end

lambda_plus = [(e_val(1)+abs(e_val(1)))*0.5 0 0; ...
          0 (e_val(2)+abs(e_val(2)))*0.5 0 ; ...
          0 0 (e_val(3)+abs(e_val(3)))*0.5];
     

lambda_minus = [(e_val(1)-abs(e_val(1)))*.5 0 0; ...
          0 (e_val(2)-abs(e_val(2)))*0.5 0 ; ...
          0 0 (e_val(3)-abs(e_val(3)))*0.5];      
      
A_plus = S_inv * Ca_inv * lambda_plus * Ca * S;      
A_minus = S_inv * Ca_inv * lambda_minus * Ca * S;   


A=A_plus - A_minus;

end

function [vhr, vhl] =boundariesWeno (p, inp)
% Interpolates the function p with weno algorithm at i+1/2 and i-1/2 points

vhr=zeros(inp.SIZEX,1);
vhl=vhr;

vhr_dummy=zeros(inp.N_N_POINTS+1,1);
vhl_dummy=vhr_dummy;

for k=inp.V1:inp.VN
    f(:,1)=p(k-inp.N_N_POINTS : k+inp.N_N_POINTS);
    if isreal(f)==0
       error ('imaginario') 
    end
    
    for l=1:inp.TOTAL_POINTS
        
        n = l+inp.N_N_POINTS;        
        
        %Calculo de v- +1/2
        xhr = inp.TOTAL_POINTS - l + 1.5;
        vhr_dummy(l) = lagrange( f(l:n), inp.X_INTERP, inp.TOTAL_POINTS, xhr);
        
        %Calculo de v+ -1/2        
        xhl = inp.TOTAL_POINTS - l + 0.5;
        vhl_dummy(l) = lagrange( f(l:n), inp.X_INTERP, inp.TOTAL_POINTS, xhl);
                    
    end
    
    beta=smoothFactor(f, inp.TOTAL_POINTS);
    C = [1/16 5/8 5/16];
    C_l = [5/16 5/8 1/16];
    
    alpha=C'./((inp.EPSILON+beta).^2);
    alpha_l = C_l' ./ (inp.EPSILON+beta).^2;
    
    w=alpha./sum(alpha); 
    w_l = alpha_l./sum(alpha_l);
    
    for m=1:inp.TOTAL_POINTS
        vhr(k)=vhr(k)+w(m)*vhr_dummy(m);
        vhl(k)=vhl(k)+w_l(m)*vhl_dummy(m);
        if w(m) < 0 
           error ('negative weight') 
        end
    end
    
end

% Boundary conditions
for i=1:inp.V1-1  
    
    vhr(i)=p(1);
    vhr(i+inp.VN)=p(inp.VN);
   
    vhl(i)=p(1);
    vhl(i+inp.VN)=p(inp.VN);

end



end

function [beta]=smoothFactor(f, s)
% Calculates smooth factor for third or forth order weno interpolation

    if s == 2
        beta=[ (f(2)-f(1))^2 ; (f(3)-f(2))^2 ];
    end
        
    if s == 3
        beta = [ 1/3*( 4*f(1)^2 - 19*f(1)*f(2) + 25*f(2)^2 + 11*f(1)*f(3) - 31*f(2)*f(3) +10*f(3)^2) ;...
            1/3*( 4*f(2)^2 -13*f(2)*f(3) + 13*f(3)^2 + 5*f(2)*f(4) - 13*f(3)*f(4) +4*f(4)^2 ) ;...
            1/3*( 10*f(3)^2 -31*f(3)*f(4) + 25*f(4)^2 + 11*f(3)*f(5) - 19*f(4)*f(5) +4*f(5)^2) ];
 
    end

end

function [result] = lagrange(f, x, k, n)

% Calculates a lagrange poynomial
% f = y axis values
% x = x axis values
% k = size of f and x
% n = point where to interpolate 
% inputs - f, k, x e c

result=0;
l=ones(1,k+1);

for j=1:k
    
    for i=1:k
       if i ~= j
          
           l(j) = l(j) * ( n - x(i) ) / ( x(j) - x(i) );
       end
        
    end
    result = result + f(j)*l(j);
    
end

end
