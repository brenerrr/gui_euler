function [P_final] = P_NewtonRaphson (gama, Pl, Pr, cl, cr, Po, max_i)

P=Po;
i=1;
alpha=(gama+1)/(gama-1);

while i <= max_i
    
    if i == max_i
        error ('No convergence')
    end
    
    f = (P-1).*cr./((1+alpha.*P).^0.5.*sqrt(gama*(gama-1)/2)) - 2.*cl./(gama-1).*(1-(P.*Pr./Pl).^((gama-1)/(2*gama)));
    fb = ((P-0.001)-1).*cr./((1+alpha.*(P-0.001)).^0.5.*sqrt(gama*(gama-1)/2)) - 2.*cl./(gama-1).*(1-((P-0.001).*Pr./Pl).^((gama-1)/(2*gama)));

    df=(f-fb)/(0.001);
    
    P_next = P - f/df;
    
    if abs(P_next-P)<0.0001*P
        i = max_i+1;
        P_final = P_next;
    end
    
    P = P_next;
    i=i+1;
    
  
end
