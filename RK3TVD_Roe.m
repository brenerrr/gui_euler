function [U_final, F_final]=RK3TVD_Roe(inp, U, F)
% Third Order TVD Runge Kutta 

U1=U; % Gets contour conditions
U2=U;
U_final=U;
F_final=F;% Gets contour conditions

dU_dt=-Roe(inp, U, F);

%% First substep 

U1(:,inp.V1:inp.VN)=U(:,inp.V1:inp.VN) + inp.DELTA_T*dU_dt;

aux(1,:) = U1(2,:) ./ U1(1,:); % u 
aux(2,:) = (inp.GAMA - 1) * (U1(3,:) - 0.5*U1(1,:).*aux(1,:).^2); % P

F_final(1,inp.V1:inp.VN) = U1(2,inp.V1:inp.VN); % Uses F_final to save memory
F_final(2,inp.V1:inp.VN) = U1(2,inp.V1:inp.VN) .* aux(1,inp.V1:inp.VN) + aux(2,inp.V1:inp.VN);
F_final(3,inp.V1:inp.VN) = ( U1(3,inp.V1:inp.VN) + aux(2,inp.V1:inp.VN) ) .* aux(1, inp.V1:inp.VN);

dU_dt=-Roe(inp, U1, F_final);

%% Second substep

U2(:,inp.V1:inp.VN)= 3/4 * U(:,inp.V1:inp.VN) + 1/4*U1(:,inp.V1:inp.VN) + 1/4*inp.DELTA_T*dU_dt;

aux(1,:) = U2(2,:) ./ U2(1,:); % u 
aux(2,:) = (inp.GAMA - 1) * (U2(3,:) - 0.5*U2(1,:).*aux(1,:).^2); % P

F_final(1,inp.V1:inp.VN) = U2(2,inp.V1:inp.VN);
F_final(2,inp.V1:inp.VN) = U2(2,inp.V1:inp.VN) .* aux(1,inp.V1:inp.VN) + aux(2,inp.V1:inp.VN);
F_final(3,inp.V1:inp.VN) = ( U2(3,inp.V1:inp.VN) + aux(2,inp.V1:inp.VN) ) .* aux(1, inp.V1:inp.VN);

dU_dt=-Roe(inp, U2, F_final);

%% Final substep

U_final(:,inp.V1:inp.VN)= 1/3 * U(:,inp.V1:inp.VN) + 2/3 * U2(:,inp.V1:inp.VN) + 2/3 * inp.DELTA_T * dU_dt;

aux(1,:) = U_final(2,:) ./ U_final(1,:); % u
aux(2,:) = (inp.GAMA - 1) * (U_final(3,:) - 0.5*U_final(1,:).*aux(1,:).^2); % P

F_final(1,inp.V1:inp.VN) = U_final(2,inp.V1:inp.VN);
F_final(2,inp.V1:inp.VN) = U_final(2,inp.V1:inp.VN) .* aux(1,inp.V1:inp.VN) + aux(2,inp.V1:inp.VN);
F_final(3,inp.V1:inp.VN) = ( U_final(3,inp.V1:inp.VN) + aux(2,inp.V1:inp.VN) ) .* aux(1, inp.V1:inp.VN);

end