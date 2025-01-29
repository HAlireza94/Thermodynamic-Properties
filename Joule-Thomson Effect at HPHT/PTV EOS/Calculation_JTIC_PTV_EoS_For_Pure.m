function [Pr, Tr]=Calculation_JTIC_PTV_EoS_For_Pure(p,Pc,Tc,T,w,Zc)
  R=8.3144598;


     Oa=(0.66121-(0.76105*Zc));
 Ob=(0.02207+(0.20868*Zc));
   Oc=(0.57765-(1.87080*Zc));

    a=(((Oa*(R^2)*(Tc^2)))/Pc);
   b=((Ob*R*Tc)/Pc);
 c=(Oc*R*Tc)/Pc;
    k=0.46283+(3.58230*w*Zc)+(8.19417*(w^2)*(Zc^2));       

    x0=(R*(T-(0.1*Tc)))/p;
    Tr=T/Tc;
    
    
        

        alpha=((1+(k*(1-(Tr^0.5))))^2);
        first_order_alpha=-(k/T)*(sqrt(Tr))*(sqrt(alpha));

    
    
    
      F=@(x) (T*((R/(x-b))-((a*first_order_alpha)/(((x*(x+b))+(c*(x-b)))))))...
        +(x*((((a*alpha*(b+c+(2*x)))/(((((x*(x+b))+(c*(x-b)))))^2)))-...
        ((R*T)/((x-b)^2))));
    V=fzero(F,x0);
    
    P=((R*T)/(V-b))-((a*alpha)/((((V*(V+b))+(c*(V-b))))));
    
    Pr=P/Pc;
    
 
end        