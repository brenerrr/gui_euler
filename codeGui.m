
%% Pre processing

% Ratios for sod's shock tube solution
Pr = 0.1 ; rhor = .125;
Pl = Pr*ratio; rhol = Pl;

analitica=1;

% Spatial and temporal parameters for different inicial conditions
if inp.IC==1
    INFX = -5; 
    SUPX = 5;
    INFT=0;
    SUPT=2;
    
    % Other artificial viscosity parameters that were not fully implemented
    inp.BULK=1*0;
    inp.cMASS=1e-3*0.5*0;
    
else
    INFX=0;
    SUPX=1;
    INFT=0; 
    if inp.IC~=2
        SUPT=0.3;
    else
        SUPT=0.178;
    end
    
    inp.BULK=1*0;
    inp.cMASS=0;
end

% WENO parameters 
if inp.WENO==1
    inp.TOTAL_POINTS=3;
    inp.N_N_POINTS=inp.TOTAL_POINTS-1; 
    inp.TABLE=load('ddcoefficients_withfactorial.txt'); 
    inp.C=coefficients_weno(inp.TOTAL_POINTS);  
    inp.EPSILON=1e-8;   
    inp.X_INTERP = 1:inp.N_N_POINTS*2+1; 
    for r = 1 : inp.TOTAL_POINTS+1
        for j = 1 : inp.TOTAL_POINTS
            inp.CRJ(r,j) = crj(r-2,j-1,inp.TOTAL_POINTS);
        end
    end
end

% Other parameters 
inp.CCL=1; 
inp.CCR=1;
ALPHA=0.2;
BETA=0;
D=0;
inp.GAMA=1.4; 
inp.R=287.0530; 
inp.gama=1.4;
inp.EC=8;

%% Grids

% Spatial grid
inp.SIZEX_O=round((SUPX-INFX)/inp.DELTA_X)+1;
inp.SIZEX=inp.SIZEX_O+inp.EC;
x=linspace(INFX-inp.EC/2*inp.DELTA_X,SUPX+inp.EC/2*inp.DELTA_X,inp.SIZEX);
inp.V1=inp.EC/2+1;
inp.VN=inp.SIZEX-inp.EC/2;

% Temporal grid
if inp.IC==1
    Co = (inp.GAMA*Pl./rhol)^0.5;
    inp.DELTA_T = CFL*inp.DELTA_X / ( 0 + Co );
else
    Co = (inp.GAMA*10.3333./3.857143)^0.5;

    inp.DELTA_T = CFL*inp.DELTA_X / ( 2.629369 + Co );
end
inp.SIZET=round((SUPT-INFT)/inp.DELTA_T)+1;
t=linspace(INFT,SUPT,inp.SIZET);

% Initialization of variables

U=zeros(3,inp.SIZEX,inp.SIZET); % rho / rho v / e 

F=U; % rho v / rho v² + P / (e+P)v


%% Initial Condition
Po=zeros(1,inp.SIZEX);

% Sod's shock tube
if inp.IC == 1
    U(1,1:round(inp.SIZEX/2),1)=rhol;
    U(1,round(inp.SIZEX/2):end,1)=rhor;
    
    
    Po(1:round(inp.SIZEX/2))=Pl;
    Po(round(inp.SIZEX/2):inp.SIZEX)=Pr;
    To = Po ./ (U(1,:,1) * inp.R);
    
    U(3,:,1)=Po ./ (inp.gama-1);
    
    F(2,1:round(inp.SIZEX/2),1)=Po(1:round(inp.SIZEX/2));
    F(2,round(inp.SIZEX/2):end,1)=Po(1,round(inp.SIZEX/2):end,1);
end

% Shu Osher Problem
if inp.IC == 2
    lambda = 1/8;
    
    i=0;
    flag=0;
    while flag~=1
        i=i+1;
        if x(i)<=1/8 && x(i+1)>=1/8
            div=i;
            flag=1;
        end
    end
    
    U(1,1:div-1,1) = 3.857143;
    U(1,div:end,1) = 1+0.2*sin(1/lambda*2*pi*x(div:end));
    
    veloci(1:div-1) = 2.629369;
    veloci(div:inp.SIZEX) = 0;
    
    U(2,:,1) = U(1,:,1).*veloci;
    
    Po(1:div-1) = 10.3333;
    Po(div:end) = 1;
    
    U(3,:,1) = Po ./ (inp.GAMA-1) + 0.5*U(2,:,1).*veloci;
    
    F(1,:,1) = U(2,:,1);
    
    F(2,:,1) = U(2,:,1).*veloci + Po;
    
    F(3,:,1) = (U(3,:,1) + Po).* veloci;
end

% Reflected Shu Osher Problem
if inp.IC==3 
    veloci = zeros(1,inp.SIZEX);
    Po=veloci;
    lambda = 1/8;
    
    i=0;
    flag=0;
    while flag~=1
        i=i+1;
        if x(i)<=1/8 && x(i+1)>=1/8
            div=i;
            flag=1;
        end
    end
    
    U(1,1:div-1,1) = 3.857143;
    U(1,inp.SIZEX-div:inp.SIZEX,1) = 3.857143/2*2;
    U(1,div:inp.SIZEX-div-1,1) = 1+0.2*sin(1/lambda*2*pi*x(div:inp.SIZEX-div-1));
    
    veloci(:) = 0;
    veloci(1:div-1) = 2.629369;
    veloci(inp.SIZEX-div:inp.SIZEX) = -2.629369;
        
    U(2,:,1) = U(1,:,1).*veloci;
    
    Po(:) = 1;
    Po(1:div-1) = 10.3333;
    Po(inp.SIZEX-div:inp.SIZEX) = 10.3333/2*2;    
        
    U(3,:,1) = Po ./ (inp.GAMA-1) + 0.5*U(2,:,1).*veloci;
    
    F(1,:,1) = U(2,:,1);
    F(2,:,1) = U(2,:,1).*veloci + Po;
    F(3,:,1) = (U(3,:,1) + Po).* veloci;
    
end

% Partially reflected Shu Osher Problem
if inp.IC==4
    veloci = zeros(1,inp.SIZEX);
    Po=veloci;
    lambda = 1/8;
    
    i=0;
    flag=0;
    while flag~=1
        i=i+1;
        if x(i)<=1/8 && x(i+1)>=1/8
            div=i;
            flag=1;
        end
    end
    
    U(1,1:div-1,1) = 3.857143;
    U(1,inp.SIZEX-div:inp.SIZEX,1) = 3.857143/2;
    U(1,div:inp.SIZEX-div-1,1) = 1+0.2*sin(1/lambda*2*pi*x(div:inp.SIZEX-div-1));
    
    veloci(:) = 0;
    veloci(1:div-1) = 2.629369;
    veloci(inp.SIZEX-div:inp.SIZEX) = -2.629369;
        
    U(2,:,1) = U(1,:,1).*veloci;
    
    Po(:) = 1;
    Po(1:div-1) = 10.3333;
    Po(inp.SIZEX-div:inp.SIZEX) = 10.3333/2;    
        
    U(3,:,1) = Po ./ (inp.GAMA-1) + 0.5*U(2,:,1).*veloci;
    
    F(1,:,1) = U(2,:,1);
    F(2,:,1) = U(2,:,1).*veloci + Po;
    F(3,:,1) = (U(3,:,1) + Po).* veloci;
    
end

%% Boundary points remain constant  
for n=1:inp.SIZET
    U(:,:,n)=U(:,:,1);
    F(:,:,n)=F(:,:,1);
end

%% Calculation

flagOk=1;
h= waitbar(0,'Calculating ...');
for n=1:inp.SIZET-1

    waitbar(n/(inp.SIZET-1),h,'Calculating ...');
    if fComp==1 && fRoe==0
        [U(:,:,n+1), F(:,:,n+1)]=RK3TVD_Compacto(inp, U(:,:,n), F(:,:,n));
    end
    
    if fRoe==1 && fComp==0
        [U(:,:,n+1), F(:,:,n+1)]=RK3TVD_Roe(inp, U(:,:,n), F(:,:,n));
    end
    
    if fRoe==0 && fComp==0 && inp.WENO==1
        [U(:,:,n+1), F(:,:,n+1)]=RK3TVD_WENO(inp, U(:,:,n), F(:,:,n));
    end
    
    U(:,inp.VN-3:inp.SIZEX,n+1) = U(:,inp.VN-3:inp.SIZEX,1);
    
    % Failsafe in case of NaN
    if isnan(U(1,20,n))==1
        flagOk = 0;
        n = inp.SIZET+1;
    end
end

% Compute Analytical Solution of Sod's Shock Tube
if handles.IC==1 && flagDone == 0
   anaSod (rhol,Pl,rhor,Pr,inp.DELTA_T, SUPT, SUPX, INFX) 
   flagDone=1;
end

close(h)
%% Post Processing

% Save results
if flagOk==1
    
    if fComp == 1 && fRoe ==0 && inp.WENO ==0
        name{1} = 'Com';
    end
    
    if fComp == 0 && fRoe == 1 && inp.WENO == 0
        name{1} = 'Roe';
        inp.SHEAR=0;
    end
    
    if fComp == 0 && fRoe == 0 && inp.WENO == 1
        name{1} = 'Wen';
        inp.SHEAR=0;
    end
    
    if fComp ==0 && fRoe == 1 && inp.WENO == 1
        name{1} = 'Rwe';
        inp.SHEAR=0;
    end
    
    name{2} = { 'Sod', 'Shu', 'Rsh', 'Prs' };
    name{2} = name{2}(inp.IC);
    
    name{3} =  'Cfl' ;
    
    name{4} =  'Dx' ;
    
    name{5} =  'Vis' ;
    
    name{6} =  'Rat' ;
    
    
    filename = strcat(name(1), ...
        name{2}(1),...
        name(3), num2str(CFL), ...
        name(4), num2str(inp.DELTA_X), ...
        name(5), num2str(inp.SHEAR), ...
        name(6), num2str(ratio));
    
    filename = strrep(filename, '.', '');
    
    
    save(strcat(filename{1},'.mat'), 'U', 'x', '-v7.3') ;
    
end

