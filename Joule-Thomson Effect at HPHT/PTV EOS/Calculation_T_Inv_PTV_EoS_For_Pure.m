function T_inv=Calculation_T_Inv_PTV_EoS_For_Pure(w,Tc,Pc,T_B,Zc)

  R=8.3144598;
    Oa=(0.66121-(0.76105*Zc));
 Ob=(0.02207+(0.20868*Zc));
%    Oc=(0.57765-(1.87080*Zc));

    a=(((Oa*(R^2)*(Tc^2)))/Pc);
   b=((Ob*R*Tc)/Pc);
%  c=(Oc*R*Tc)/Pc;
    k=0.46283+(3.58230*w*Zc)+(8.19417*(w^2)*(Zc^2));       


       
       F=@(x)  -(((((a*(-(k/x)*(sqrt(x/Tc))*...
                (sqrt(((1+(k*(1-((x/Tc)^0.5))))^2)))))*R*x)-...
                (R*a*((1+(k*(1-((x/Tc)^0.5))))^2))))/...
                ((R*x)^2))-((b-((a*((1+(k*(1-((x/Tc)^0.5))))^2))/(R*x)))/x);

       T_inv=fzero(F,T_B);

        

  
end