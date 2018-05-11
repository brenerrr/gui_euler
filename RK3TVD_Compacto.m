function [U_final, F_final]=RK3TVD_Compacto(inp, U, F)
% Third Order TVD Runge Kutta for 10th order compact scheme

U1=U; % Get boundaries
U2=U;
U_final=U;
F_final=F;


dU_dt(1,:)=-compactDfTenthOrder(F(1,:),inp);
dU_dt(2,:)=-compactDfTenthOrder(F(2,:),inp);
dU_dt(3,:)=-compactDfTenthOrder(F(3,:),inp);

%% First substep
U1(:,inp.V1:inp.VN)=U(:,inp.V1:inp.VN) + 1*inp.DELTA_T*dU_dt;

U1(1,inp.V1:inp.VN)=compactFilter(U1(1,inp.V1:inp.VN), inp.SIZEX_O);
U1(2,inp.V1:inp.VN)=compactFilter(U1(2,inp.V1:inp.VN), inp.SIZEX_O);
U1(3,inp.V1:inp.VN)=compactFilter(U1(3,inp.V1:inp.VN), inp.SIZEX_O);

aux(1,:) = U1(2,:) ./ U1(1,:); % u 
aux(2,:) = (inp.GAMA - 1) * (U1(3,:) - 0.5*U1(1,:).*aux(1,:).^2); % P

F_final(1,inp.V1:inp.VN) = U1(2,inp.V1:inp.VN);
F_final(2,inp.V1:inp.VN) = U1(2,inp.V1:inp.VN) .* aux(1,inp.V1:inp.VN) + aux(2,inp.V1:inp.VN);
F_final(3,inp.V1:inp.VN) = ( U1(3,inp.V1:inp.VN) + aux(2,inp.V1:inp.VN) ) .* aux(1, inp.V1:inp.VN);

dU_dt(1,:)=-compactDfTenthOrder(F_final(1,:),inp);
dU_dt(2,:)=-compactDfTenthOrder(F_final(2,:),inp);
dU_dt(3,:)=-compactDfTenthOrder(F_final(3,:),inp);
%% Second substep

U2(:,inp.V1:inp.VN)= 3/4 * U(:,inp.V1:inp.VN) + 1/4*U1(:,inp.V1:inp.VN) + 1/4*inp.DELTA_T*dU_dt;

U2(1,inp.V1:inp.VN)=compactFilter(U2(1,inp.V1:inp.VN), inp.SIZEX_O);
U2(2,inp.V1:inp.VN)=compactFilter(U2(2,inp.V1:inp.VN), inp.SIZEX_O);
U2(3,inp.V1:inp.VN)=compactFilter(U2(3,inp.V1:inp.VN), inp.SIZEX_O);

aux(1,:) = U2(2,:) ./ U2(1,:); % u 
aux(2,:) = (inp.GAMA - 1) * (U2(3,:) - 0.5*U2(1,:).*aux(1,:).^2); % P

F_final(1,inp.V1:inp.VN) = U2(2,inp.V1:inp.VN);
F_final(2,inp.V1:inp.VN) = U2(2,inp.V1:inp.VN) .* aux(1,inp.V1:inp.VN) + aux(2,inp.V1:inp.VN);
F_final(3,inp.V1:inp.VN) = ( U2(3,inp.V1:inp.VN) + aux(2,inp.V1:inp.VN) ) .* aux(1, inp.V1:inp.VN);

dU_dt(1,:)=-compactDfTenthOrder(F_final(1,:),inp);
dU_dt(2,:)=-compactDfTenthOrder(F_final(2,:),inp);
dU_dt(3,:)=-compactDfTenthOrder(F_final(3,:),inp);

%% Final substep

U_final(:,inp.V1:inp.VN)= 1/3 * U(:,inp.V1:inp.VN) + 2/3 * U2(:,inp.V1:inp.VN) + 2/3 * inp.DELTA_T * dU_dt;

U_final(1,inp.V1:inp.VN)=compactFilter(U_final(1,inp.V1:inp.VN), inp.SIZEX_O);
U_final(2,inp.V1:inp.VN)=compactFilter(U_final(2,inp.V1:inp.VN), inp.SIZEX_O);
U_final(3,inp.V1:inp.VN)=compactFilter(U_final(3,inp.V1:inp.VN), inp.SIZEX_O);

aux(1,:) = U_final(2,:) ./ U_final(1,:); % u
aux(2,:) = (inp.GAMA - 1) * (U_final(3,:) - 0.5*U_final(1,:).*aux(1,:).^2); % P

F_final(1,inp.V1:inp.VN) = U_final(2,inp.V1:inp.VN);
F_final(2,inp.V1:inp.VN) = U_final(2,inp.V1:inp.VN) .* aux(1,inp.V1:inp.VN) + aux(2,inp.V1:inp.VN);
F_final(3,inp.V1:inp.VN) = ( U_final(3,inp.V1:inp.VN) + aux(2,inp.V1:inp.VN) ) .* aux(1, inp.V1:inp.VN);

[F_final(1,inp.V1:inp.VN), F_final(2,inp.V1:inp.VN), F_final(3,inp.V1:inp.VN)] = a_viscosity(U_final, F_final, inp);
end

function [F_new1, F_new2, F_new3] = a_viscosity(U, F, inp)
% Adds an artificial viscosity into the fluxes of linear momentum and
% energy. Currently it is also implemented a numerical diffusivity into the
% density flux, but this is not validated yet, what makes it unreliable.
    
    u = U(2,:) ./ U(1,:); % Velocity
    P = (inp.GAMA - 1) * (U(3,:) - 0.5*U(2,:).*u); % Pressure
    
    % Not validated yet ***************************************************
    cp = 1;
    T = P./(U(1,:)*inp.R);
    dT = (T(inp.V1+1:inp.VN+1)-T(inp.V1-1:inp.VN-1))/inp.DELTA_X/2;
    dP = (P(inp.V1+1:inp.VN+1)-P(inp.V1-1:inp.VN-1))/inp.DELTA_X/2;
    ds = zeros(1, inp.SIZEX);
    ds(inp.V1:inp.VN) = cp./T(inp.V1:inp.VN) .* dT + inp.R./P(inp.V1:inp.VN) .* dP;
    ds(1:inp.V1-1)=ds(inp.V1);
    ds(inp.VN+1:inp.SIZEX)=ds(inp.VN); 
    ds = explicit_1deriv_fourthorder(ds, inp);

    % Used for contact discontinuities  and temperature gradients
    chi = inp.cMASS * mean(sqrt(inp.GAMA*P(inp.V1:inp.VN)./U(1,(inp.V1:inp.VN))))/cp *inp.DELTA_X^6 * ...
        gaussian_filter(abs(ds), inp.SIZEX_O);
    
    chi = chi .* (U(1,inp.V1+1:inp.VN+1)-U(1,inp.V1-1:inp.VN-1))/inp.DELTA_X/2;
    % Not validated yet ***************************************************   
    
    du_dx = (u(inp.V1+1:inp.VN+1)-u(inp.V1-1:inp.VN-1))/inp.DELTA_X/2;

    f=zeros(inp.SIZEX,1);

    f(inp.V1:inp.VN) = dvfSecondOrder (u, inp); % Explicit Fith derivative
    f(1:inp.V1-1)=f(inp.V1);
    f(inp.VN+1:inp.SIZEX)=f(inp.VN);    
    f=U(1,:).*gaussian_filter(abs(f'),inp.SIZEX);

    nr = inp.DELTA_X^6.*f(inp.V1:inp.VN);

    shear = inp.SHEAR*nr;
    bulk = inp.BULK*nr;

    tau = shear*2.*du_dx + (bulk - 2/3*shear).*du_dx;

    F_new1 = F(1,inp.V1:inp.VN) - chi;
    F_new2 =F(2,inp.V1:inp.VN)-tau;
    F_new3= ( U(3,inp.V1:inp.VN) + P(inp.V1:inp.VN) - tau ) .* u(inp.V1:inp.VN);
end

function [df] = dvfSecondOrder (f, inp)
% Calculates the fith derivative with second order precision

      df= ( -1/2*(f(inp.V1-3:inp.VN-3)-f(inp.V1+3:inp.VN+3)) + ...
              2* (f(inp.V1-2:inp.VN-2)-f(inp.V1+2:inp.VN+2)) - ...
              5/2*(f(inp.V1-1:inp.VN-1)-f(inp.V1+1:inp.VN+1)) + ...
              0*f(inp.V1:inp.VN) ) / inp.DELTA_X^5;
end

function [f_filtered] = gaussian_filter(f, sizef)
% Applies a gaussian filter to function f

    f_filtered=f;
    j=0;    
    for i=5:sizef-4   
       j=j+1;
       f_filtered(i) = 3565/10368*f(i) + 3091/12960*(f(i-1)+f(i+1)) + 1997/25920*(f(i-2)+f(i+2)) + ...
           149/12960*(f(i-3)+f(i+3))+107/103680*(f(i-4)+f(i+4));
    end
    
end

function [f_filtered] = compactFilter(f, sizef)
% Smooths function f with a compact filter that dumps only very high wave
% numbers
    
    alpha=0.66624;
    beta=0.16688;
    a = 0.99965;
    b = 0.66652;
    c = 0.16674;
    d = 4e-5;
    e = -5e-6;
    sizeb=sizef-8;
    B=zeros(sizeb,1);
    
    j=0;
    for i=5:sizef-4
        j=j+1;
        B(j) = a*f(i) + b*( f(i-1)+f(i+1) ) +c*( f(i-2)+f(i+2) ) +d*( f(i-3)+f(i+3) ) +...
            e*( f(i-4)+f(i+4));
    end
    
    A = pentadiagonal (beta,alpha,1,alpha,beta,sizef-8);
    
    f_filtered=f;
    
    B(1)=B(1) - (beta*f_filtered(3) + alpha*f_filtered(4));
    B(2)=B(2) - beta*f_filtered(4);
    B(sizeb) = B(sizeb) - (alpha*f_filtered(sizef-3) + beta*f_filtered(sizef-2));
    B(sizeb-1) = B(sizeb-1) - beta*f_filtered(sizef-3);

    
    f_filtered(5:sizef-4)=linsolve(A,B);

end

function [df] = compactDfTenthOrder (f,inp)
% Calculates the derivative of f with tenth order of precision with
% constant boundaries conditions (f' at boundaries = 0)

    alpha=0.5;
    beta=1/20;
    a=17/12;
    b=101/150;
    c=1/100;
    
    B = c/(6*inp.DELTA_X)*(f(inp.V1+3:inp.VN+3)-f(inp.V1-3:inp.VN-3)) + ...
        b/(4*inp.DELTA_X)*(f(inp.V1+2:inp.VN+2)-f(inp.V1-2:inp.VN-2)) + ...
        a/(2*inp.DELTA_X)*(f(inp.V1+1:inp.VN+1)-f(inp.V1-1:inp.VN-1));
    
    
    A=pentadiagonal(beta,alpha,1,alpha,beta,inp.SIZEX_O);
    
    df=linsolve(A,B');
end

