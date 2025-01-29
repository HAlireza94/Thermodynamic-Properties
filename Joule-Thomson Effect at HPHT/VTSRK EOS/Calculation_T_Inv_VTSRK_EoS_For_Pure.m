function T_inv=Calculation_T_Inv_VTSRK_EoS_For_Pure(alpha_Function,w,Tc,Pc,T_B,Zc)

  R=8.3144598;
  a=(((0.42747*(R^2)*(Tc^2)))/Pc);
  b=(0.08664*R*Tc)/Pc;
c1=(-45.7247*((1/3)-Zc));
  c2=((-2.184*exp(c1))+0.2658);
  c=((((1/3)-Zc)*((R*Tc)/Pc))*c2);


  switch lower(alpha_Function)
    
    case 'original'
        k=0.480+(1.574*w)-(0.176*w^2);

        F=@(x)  -(((((a*(-(k/x)*(sqrt(x/Tc))*(sqrt(((1+(k*(1-((x/Tc)^0.5))))^2)))))*R*x)-...
                (R*a*((1+(k*(1-((x/Tc)^0.5))))^2))))/...
    ((R*x)^2))-((b-c-((a*((1+(k*(1-((x/Tc)^0.5))))^2))/(R*x)))/x);

     T_inv=fzero(F,T_B);
                
    case 'ozokwelu'
          
         m=0.266+(0.4459*w^0.5);
         n2=(1/m)*(0.2469+(0.7495*w));
        F=@(x)  -(((((a*(((-2*m*n2)/Tc)*((x/Tc)^(n2-1))*(exp(2*m*(1-((x/Tc)^n2))))))*R*x)-...
        (R*a*(((exp(m*(1-((x/Tc)^n2))))^2)))))/...
        ((R*x)^2))-((b-c-((a*((exp(m*(1-((x/Tc)^n2))))^2))/(R*x)))/x);

        T_inv=fzero(F,T_B);

        
           case 'trebble'
         
         if w<-0.1
             
             q1=0.66208+(4.63961*w)+(7.45183*w^2);
             
         elseif w>=-0.1 && w<=0.4
             
             q1=0.35+(0.7924*w)+(0.1875*w^2)-(28.93*((0.3-Zc)^2));
             
         elseif w>0.4
             
             q1=0.32+(0.9424*w)-(28.93*((0.3-Zc)^2));
             
         end
        
        
       F=@(x)  -(((((a*(-(q1/Tc)*(exp(q1*(1-(x/Tc))))))*R*x)-...
                (R*a*(exp(q1*(1-(x/Tc)))))))/...
    ((R*x)^2))-((b-c-((a*(exp(q1*(1-(x/Tc)))))/(R*x)))/x);

     T_inv=fzero(F,T_B);

        
          case 'nasrifar'
         
%                  if Tr>1
                    
                k=0.480+(1.574*w)-(0.175*w^2);
                b1=0.25*(12-(11*k)+(k^2));
                b2=0.5*(-6+(9*k)-(k^2));
                b3=0.25*(4-(7*k)+(k^2));                                 

        F=@(x)  -(((((a*(-((Tc*b1)/(x^2))-...
        ((2*(Tc^2)*b2)/(x^3))-...
        ((3*(Tc^3)*b3)/(x^4))))*R*x)-...
                (R*a*((((b1/(x/Tc))+(b2/...
        ((x/Tc))^2)+(b3/((x/Tc))^3)))))))/...
    ((R*x)^2))-((b-c-((a*((((b1/(x/Tc))+(b2/...
        ((x/Tc))^2)+(b3/((x/Tc))^3)))))/(R*x)))/x);

     T_inv=fzero(F,T_B);
   
%                 elseif Tr<1 || Tr==1                      
%             
%          
%            
%         k=0.480+(1.574*w)-(0.176*w^2);
%         alpha=((1+(k*(1-(Tr^0.5))))^2);
%         first_order_alpha=-(k/T)*(sqrt(Tr))*(sqrt(alpha));
%                         
%                 end

     
     
  end

  
end