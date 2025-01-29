function [Pr, Tr]=PR_Equation_of_Satate(alpha_Function,p,Pc,Tc,T,w)
R=83.14472;

  a=(0.45724*(R^2)*(Tc^2))/Pc;
  b=(0.07780*R*Tc)/Pc;

    x0=(R*T)/p;
    Tr=T/Tc;
    



switch lower(alpha_Function)
    
    case 'coquelet'
        
    k=0.41287+(1.34494*w)+(0.00421*w^2);

        alpha=exp(k*(1-Tr));
    first_order_alpha=-(k/Tc)*alpha;
    
    F1=@(x) (T*((R/(x-b))-((a*first_order_alpha)/((x^2)+(2*b*x)-(b^2)))))...
        +(x*((((2*a*alpha*(x+b))/(((x^2)+(2*b*x)-(b^2))^2)))-...
        ((R*T)/((x-b)^2))));
    
    V=fzero(F1,x0);
    
    P=((R*T)/(V-b))-((a*alpha)/((V^2)+(2*b*V)-b^2));
    
    Pr=P/Pc;
    
    case 'original'
        
     k=0.37464+(1.54226*w)-(0.26992*w^2);

        alpha=((1+(k*(1-(Tr^0.5))))^2);
    first_order_alpha=-(k/T)*(sqrt(Tr))*(sqrt(alpha));
    
    F1=@(x) (T*((R/(x-b))-((a*first_order_alpha)/((x^2)+(2*b*x)-(b^2)))))...
        +(x*((((2*a*alpha*(x+b))/(((x^2)+(2*b*x)-(b^2))^2)))-...
        ((R*T)/((x-b)^2))));
    
    V=fzero(F1,x0);
    
    P=((R*T)/(V-b))-((a*alpha)/((V^2)+(2*b*V)-b^2));
    
    Pr=P/Pc;
    
    case 'haghtalab'
        
    case 'saffari'
        
end

end        