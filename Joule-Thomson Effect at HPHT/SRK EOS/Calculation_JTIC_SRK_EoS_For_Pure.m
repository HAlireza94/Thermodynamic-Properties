function [Pr, Tr]=Calculation_JTIC_SRK_EoS_For_Pure(alpha_Function,p,Pc,Tc,T,w,Zc)
  
R=8.3144598;

a=(((0.42747*(R^2)*(Tc^2)))/Pc);
b=(0.08664*R*Tc)/Pc;
    x0=(R*T)/p;
    Tr=T/Tc;
    
    
    switch lower(alpha_Function)

        case 'original'
        
        k=0.480+(1.574*w)-(0.176*w^2);
        alpha=((1+(k*(1-(Tr^0.5))))^2);
        first_order_alpha=-(k/T)*(sqrt(Tr))*(sqrt(alpha));
            
        case 'ozokwelu'
             
          m=0.266+(0.4459*w^0.5);
          n2=(1/m)*(0.2469+(0.7495*w));

            alpha=((exp(m*(1-((T/Tc)^n2))))^2);
            first_order_alpha=((-2*m*n2)/Tc)*(Tr^(n2-1))*(exp(2*m*(1-(Tr^n2))));
            
            
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
     

   
                elseif Tr<1 || Tr==1                      
            
         
           
        k=0.480+(1.574*w)-(0.176*w^2);
        alpha=((1+(k*(1-(Tr^0.5))))^2);
        first_order_alpha=-(k/T)*(sqrt(Tr))*(sqrt(alpha));
                        
                end
   
        
    end
    
    
    
    F=@(x) (T*((R/(x-b))-((a*first_order_alpha)/(x*(x+b)))))...
        +(x*((((a*alpha*((2*x)+b))/((x*(x+b))^2)))-...
        ((R*T)/((x-b)^2))));
    
    V=fzero(F,x0);
    
    P=((R*T)/(V-b))-((a*alpha)/(V*(V+b)));
    
    Pr=P/Pc;
    
 
end        