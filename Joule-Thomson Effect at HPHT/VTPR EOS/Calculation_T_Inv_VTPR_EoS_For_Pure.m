function T_inv=Calculation_T_Inv_VTPR_EoS_For_Pure(alpha_Function,w,Tc,Pc,T_B,Zc)

  R=8.3144598;
  a=(0.45724*(R^2)*(Tc^2))/Pc;
  b=(0.07780*R*Tc)/Pc;
  c=(-0.252*((R*Tc)/Pc)*((1.5448*Zc)-0.4024));
  switch lower(alpha_Function)
    
    case 'coquelet'
          k=0.41287+(1.34494*w)+(0.00421*w^2);

        F=@(x)  -(((((a*(-(k/Tc)*(exp(k*(1-(x/Tc))))))*R*x)-...
                (R*a*(exp(k*(1-(x/Tc)))))))/...
    ((R*x)^2))-((b-c-((a*exp(k*(1-(x/Tc))))/(R*x)))/x);

T_inv=fzero(F,T_B);
        
    case 'original'
       k=0.37464+(1.54226*w)-(0.26992*w^2);
       
       F=@(x)  -(((((a*(-(k/x)*(sqrt(x/Tc))*...
                (sqrt(((1+(k*(1-((x/Tc)^0.5))))^2)))))*R*x)-...
                (R*a*((1+(k*(1-((x/Tc)^0.5))))^2))))/...
                ((R*x)^2))-((b-c-((a*((1+(k*(1-((x/Tc)^0.5))))^2))/(R*x)))/x);

       T_inv=fzero(F,T_B);

        
    case 'haghtalab'
          
           n2=1.7309+(1.6571*w)+(0.1554*w^2);            

        F=@(x)  -(((((a*((-(exp(1-(n2^(log((x/Tc))))))*(x^(log(n2)-1))*(log(n2)))/(Tc^(log(n2)))))*R*x)-...
        (R*a*(exp(1-(n2^(log(x/Tc))))))))/...
        ((R*x)^2))-((b-c-((a*(exp(1-(n2^(log(x/Tc))))))/(R*x)))/x);

        T_inv=fzero(F,T_B);
        
        
    case 'saffari'
        
        k1=0.003091+(0.013145*w);
            k2=-0.006478+(0.482173*w);
            k3=(3.58616*w)+0.721306;            
         F=@(x)  -(((((a*(exp((k1*(x/Tc))+(k3*(1-((x/Tc)^0.5)))))*...
                (((2*k2*(Tc^1.5)*(x^(k2-1)))-(k3*Tc*(x^(k2-0.5)))+...
                (2*k1*sqrt(Tc)*(x^k2)))/(2*(Tc^(k2+1.5)))))*R*x)-...
                   (R*a*(exp((k1*(x/Tc))+(k2*(log(x/Tc)))+(k3*(1-((x/Tc)^0.5))))))))/...
                   ((R*x)^2))-((b-c-((a*(exp((k1*(x/Tc))+(k2*(log(x/Tc)))+...
                   (k3*(1-((x/Tc)^0.5))))))/(R*x)))/x);

               
             T_inv=fzero(F,T_B);
        
        
  end

  
end