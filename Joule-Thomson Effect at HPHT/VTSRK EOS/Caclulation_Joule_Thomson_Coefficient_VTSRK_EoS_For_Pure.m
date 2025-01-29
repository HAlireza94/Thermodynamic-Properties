function [J_T, Speed_Sound, Cp, Cv]=Caclulation_Joule_Thomson_Coefficient_VTSRK_EoS_For_Pure...
    (alpha_function,Tc,Pc,w,T,P,A,B,C,D,E,Zc,a0,Mw)

% Pc=pc.*10^5;
R=8.3144598;
    Tr=T/Tc;

    x0=((R*(T+Tc+a0))/P);
    
    
    a=(((0.42747*(R^2)*(Tc^2)))/Pc);
    b=(0.08664*R*Tc)/Pc;
    c1=(-45.7247*((1/3)-Zc));
    c2=((-2.184*exp(c1))+0.2658);
    c=((((1/3)-Zc)*((R*Tc)/Pc))*c2);


 switch lower(alpha_function)

        case 'original'
        
       k=0.480+(1.574*w)-(0.176*w^2);
        alpha=((1+(k*(1-(Tr^0.5))))^2);
        first_order_alpha=-(k/T)*(sqrt(Tr))*(sqrt(alpha));
         second_order_alpha=(-k*((-k*sqrt(T)-...
                            (sqrt(Tc)*sqrt(alpha)))))/(2*Tc*(T^1.5));

         F=@(x) ((R*T)/(x+c-b))-((a*((1+(k*(1-((T/Tc)^0.5))))^2))/...
             ((x+c)*(x+c+b)))-P;
     
         V=fzero(F,x0); % OF : Volume

        case 'ozokwelu'
            
       m=0.266+(0.4459*w^0.5);
       n2=(1/m)*(0.2469+(0.7495*w));
        alpha=((exp(m*(1-((T/Tc)^n2))))^2);
        first_order_alpha=((-2*m*n2)/Tc)*(Tr^(n2-1))*(exp(2*m*(1-(Tr^n2))));
        second_order_alpha=((-2*m*n2*(T^(n2-2))*((((Tc^n2)*(n2-1))-...
                           (2*m*n2*T^n2)))*(exp(2*m*(1-(Tr^n2)))))/(Tc^(2*n2)));
        F=@(x) ((R*T)/(x+c-b))-((a*((exp(m*(1-((T/Tc)^n2))))^2))/((x+c)*(x+c+b)))-P;
     
        V=fzero(F,x0);           

        
         case 'trebble'
         
         if w<-0.1
             
             q1=0.66208+(4.63961*w)+(7.45183*w^2);
             
         elseif w>=-0.1 && w<=0.4
             
             q1=0.35+(0.7924*w)+(0.1875*w^2)-(28.93*((0.3-Zc)^2));
             
         elseif w>0.4
             
             q1=0.32+(0.9424*w)-(28.93*((0.3-Zc)^2));
             
         end

         alpha=exp(q1*(1-Tr));
         first_order_alpha=-(q1/Tc)*alpha;
        second_order_alpha=((q1/Tc)^2)*alpha;
         
         F=@(x) ((R*T)/(x+c-b))-((a*(exp(q1*(1-Tr))))/((x+c)*(x+c+b)))-P;
     
        V=fzero(F,x0);
        
         case 'nasrifar'
         
                 if Tr>1
                    
                k=0.480+(1.574*w)-(0.175*w^2);
                b1=0.25*(12-(11*k)+(k^2));
                b2=0.5*(-6+(9*k)-(k^2));
                b3=0.25*(4-(7*k)+(k^2));
                
                 alpha=(((b1/(T/Tc))+(b2/...
        ((T/Tc))^2)+(b3/((T/Tc))^3)));
    
    first_order_alpha=-((Tc*b1)/(T^2))-...
        ((2*(Tc^2)*b2)/(T^3))-...
        ((3*(Tc^3)*b3)/(T^4));
     

    second_order_alpha=((2*Tc*b1)/(T^3))+...
        ((6*(Tc^2)*b2)/(T^4))+...
        ((12*(Tc^3)*b3)/(T^5));
    
     F=@(x) ((R*T)/(x+c-b))-((a*((((b1/(T/Tc))+(b2/...
        ((T/Tc))^2)+(b3/((T/Tc))^3)))))/((x+c+b)*(x+c)))-P;
     
        V=fzero(F,x0);

                elseif Tr<1 || Tr==1                      
            
         
           
        k=0.480+(1.574*w)-(0.176*w^2);
        alpha=((1+(k*(1-(Tr^0.5))))^2);
        first_order_alpha=-(k/T)*(sqrt(Tr))*(sqrt(alpha));
        second_order_alpha=(-k*((-k*sqrt(T)-...
                            (sqrt(Tc)*sqrt(alpha)))))/(2*Tc*(T^1.5));
                        
       F=@(x) ((R*T)/(x+c-b))-((a*((1+(k*(1-(Tr^0.5))))^2))/((x+c+b)*(x+c)))-P;
     
        V=fzero(F,x0);
                        
                 end
     case 'jub'
       k=0.480+(1.5963*w)-(0.2963*w^2)+(0.1223*w^3);
        alpha=((1+(k*(1-(Tr^0.5))))^2);
        first_order_alpha=-(k/T)*(sqrt(Tr))*(sqrt(alpha));
         second_order_alpha=(-k*((-k*sqrt(T)-...
                            (sqrt(Tc)*sqrt(alpha)))))/(2*Tc*(T^1.5));

         F=@(x) ((R*T)/(x+c-b))-((a*((1+(k*(1-((T/Tc)^0.5))))^2))/...
             ((x+c)*(x+c+b)))-P;
     
         V=fzero(F,x0); % OF : Volume


 end


    Cp_star=A+(B*T)+(C*(T^2))+(D*(T^3))+(E*(T^4));
    
    Cv_star=Cp_star-R;
    
%     Tr=T/Tc;
    
%     Z=(P*OF)/(R*T);
              
      Cv=Cv_star+(T*...
     ((a*second_order_alpha)/(b))*...
     (log((b/(V+c))+1)));
    
    first_order_P_T=(R/(V+c-b))-...
        (a*first_order_alpha/((V+c)*(V+c+b)));
    
    first_order_P_V=((-R*T)/((V+c-b)^2))+...
        ((a*alpha*((2*c)+(2*V)+b)))/((((V+c)*(V+c+b)))^2);
    
    Cp=Cv-((T*(first_order_P_T^2))/(first_order_P_V));
    
     Gamma=Cp/Cv;
    
     SS=sqrt((-V^2)*(Gamma/Mw)*first_order_P_V); 
    
    Speed_Sound=SS*sqrt(1000);
    
    J_T=((T*((-first_order_P_T)/first_order_P_V))-V)/Cp;
    



end