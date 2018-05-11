function [ ] = anaSod (rhol,Pl,rhor,Pr,delta_t, supt, supx, infx)

% Got equations from Hirsch cap. 16. Expression for finding P was actually
% wrong, so I had to find one with brute force. Equation for pressure in
% expansion zone was also wrong, so I had to take it from an ALTAIR
% resolution.

% This code is messy as hell and has no comments

gama=1.4;
R=287.0530;


Ul=0;
cl=sqrt(gama*Pl/rhol);

Ur=0;
cr=sqrt(gama*Pr/rhor);

Po=2;

x=infx:0.001:supx;
size_x=length(x);
xo=(supx+infx)/2;

sizet=round(supt/delta_t)+1;

f=zeros(3,size_x,sizet);

f(1,round(1:size_x/2),1)=rhol;  
f(1,round(1:size_x/2):end,1)=rhor;
f(2,round(1:size_x/2),1)=rhol*Ul;  
f(2,round(1:size_x/2),1)=rhor*Ur;
f(3,round(1:size_x/2),1)=Pl/(gama-1);  
f(3,round(1:size_x/2):end,1)=Pr/(gama-1);


max_i=1000;
alpha=(gama+1)/(gama-1);

P=P_NewtonRaphson (gama, Pl, Pr, cl, cr, Po, max_i);

P_post=P*Pr;

rho_post=rhor * (1 + alpha*P) / (alpha + P);

U_post= 2/(gama-1)*cl*(1-(P*Pr/Pl)^((gama-1)/(2*gama)));
U_middle=U_post;

U_shock=cr^2*(P-1)/(gama*(U_post));

c_post=sqrt(gama*P_post/rho_post);

rho_middle=(P_post/Pl)^(1/gama)*rhol;

P_middle=P_post;

c_middle=sqrt(gama*P_middle/rho_middle);

for n=2:sizet  
    
    x1=xo-cl*(n-1)*delta_t;
    x2=xo+(U_post-c_middle)*(n-1)*delta_t;
    x3=xo+U_post*(n-1)*delta_t;
    x4=xo+U_shock*(n-1)*delta_t;
    
     
    for i=1:size_x
        if x(i)<x1
           f(1,i,n)=rhol;
           f(2,i,n)=rhol*Ul;
           f(3,i,n)=Pl/(gama-1)+0.5*rhol*Ul^2;      
        end
        
        if x(i)>=x1 && x(i)<x2
            U=2/(gama+1)*( (x(i)-xo)/((n-1)*delta_t) + cl );
            P=Pl*(1-(gama-1)/2*U/cl)^(2*gama/(gama-1));
            rho=rhol*(P/Pl)^(1/gama);
            
            f(1,i,n)=rho;
            f(2,i,n)=rho*U;
            f(3,i,n)=P/(gama-1)+0.5*rho*U^2;    
        end
        
        if x(i)>=x2 && x(i)<x3
            f(1,i,n)=rho_middle;
            f(2,i,n)=rho_middle*U_middle;
            f(3,i,n)=P_middle/(gama-1)+0.5*rho_middle*U_middle^2;
        end
        
        if x(i)>=x3 && x(i)<x4
            f(1,i,n)=rho_post;
            f(2,i,n)=rho_post*U_post;
            f(3,i,n)=P_post/(gama-1)+0.5*rho_post*U_post^2;
        end
        
        if x(i)>=x4
            f(1,i,n)=rhor;
            f(2,i,n)=rhor*Ur;
            f(3,i,n)=Pr/(gama-1)+0.5*rhor*Ur^2;
        end
              
    end
    
end

filename = strcat('Sod','Rat',num2str(Pl/Pr),'Reference.mat');
save(filename,'f', 'x', '-v7.3') ;

end