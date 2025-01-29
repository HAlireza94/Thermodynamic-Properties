function T_inv=calculation_T_inv(alpha_Function,w,Tc,Pc,T_B)

R=8.3144598;



  a=(0.45724*(R^2)*(Tc^2))/Pc;
  b=(0.07780*R*Tc)/Pc;

  switch lower(alpha_Function)
    
    case 'coquelet'
          k=0.41287+(1.34494*w)+(0.00421*w^2);

        F2=@(x)  -(((((a*(-(k/Tc)*(exp(k*(1-(x/Tc))))))*R*x)-...
                (R*a*(exp(k*(1-(x/Tc)))))))/...
    ((R*x)^2))-((b-((a*exp(k*(1-(x/Tc))))/(R*x)))/x);

T_inv=fzero(F2,T_B);
        
    case 'original'
             k=0.37464+(1.54226*w)-(0.26992*w^2);
 F=@(x)  -(((((a*(-(k/x)*(sqrt(x/Tc))*...
                (sqrt(((1+(k*(1-((x/Tc)^0.5))))^2)))))*R*x)-...
                (R*a*((1+(k*(1-((x/Tc)^0.5))))^2))))/...
    ((R*x)^2))-((b-((a*((1+(k*(1-((x/Tc)^0.5))))^2))/(R*x)))/x);

T_inv=fzero(F,T_B);

        
    case 'haghtalab'
        
    case 'saffari'
        
  end

  
end