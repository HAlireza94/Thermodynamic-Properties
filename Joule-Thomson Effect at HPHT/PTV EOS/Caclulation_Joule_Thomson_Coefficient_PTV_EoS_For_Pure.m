function [J_T, Speed_Sound, Cp, Cv]=...
    Caclulation_Joule_Thomson_Coefficient_PTV_EoS_For_Pure(Tc,Pc,w,T,P,A,B,C,D,E,Zc,a0,Mw)

% Pc=pc.*10^5;
R=8.3144598;
    Tr=T/Tc;

    x0=((R*(T+Tc+a0))/P);
    
    
  Oa=(0.66121-(0.76105*Zc));
 Ob=(0.02207+(0.20868*Zc));
   Oc=(0.57765-(1.87080*Zc));

    a=(((Oa*(R^2)*(Tc^2)))/Pc);
   b=((Ob*R*Tc)/Pc);
 c=(Oc*R*Tc)/Pc;
    k=0.46283+(3.58230*w*Zc)+(8.19417*(w^2)*(Zc^2));       
        
    alpha=((1+(k*(1-(Tr^0.5))))^2);
        first_order_alpha=-(k/T)*(sqrt(Tr))*(sqrt(alpha));
         second_order_alpha=(-k*((-k*sqrt(T)-...
                            (sqrt(Tc)*sqrt(alpha)))))/(2*Tc*(T^1.5));

         F=@(x)((R*T)/(x-b))-((a*((1+(k*...
             (1-((T/Tc)^0.5))))^2))/((x*(x+b))+(c*(x-b))))-P;
     
         V=fzero(F,x0); % OF : Volume
   

    Cp_star=A+(B*T)+(C*(T^2))+(D*(T^3))+(E*(T^4));
    
    Cv_star=Cp_star-R;
    
%     Tr=T/Tc;
    
%     Z=(P*OF)/(R*T);
                   SQ=sqrt((c^2)+(b^2)+(6*b*c)); %true

    Cv=Cv_star-(T*...
     ((a*second_order_alpha)/(SQ)))*...
     (log((((2*V)-SQ+c+b)/((2*V)+SQ+c+b)))); %true
    
 
 
 
    first_order_P_T=(R/(V-b))-...
        (a*first_order_alpha/((V*(V+b))+(c*(V-b)))); % true
    
    first_order_P_V=((-R*T)/((V-b)^2))+...
        ((a*alpha*(b+c+(2*V)))/...
        (((V*(V+b))+(c*(V-b)))^2)); % true
    
    Cp=Cv-((T*(first_order_P_T^2))/(first_order_P_V));
    
       Gamma=Cp/Cv;
    
     SS=sqrt((-V^2)*(Gamma/Mw)*first_order_P_V); 
    
    Speed_Sound=SS*sqrt(1000);
        
    J_T=((T*((-first_order_P_T)/first_order_P_V))-V)/Cp; % true
    



end