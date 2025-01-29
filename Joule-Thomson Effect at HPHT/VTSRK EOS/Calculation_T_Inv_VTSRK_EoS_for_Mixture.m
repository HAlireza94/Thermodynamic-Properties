function OBJ=Calculation_T_Inv_VTSRK_EoS_for_Mixture...
    (alpha_function,T1,Tc,Pc,w,Kij,Yi,N,Zc)

R=83.14472;
a=zeros(1,N);b=zeros(1,N);c1=zeros(1,N);
Tr=zeros(1,N);mb=zeros(1,N);c2=zeros(1,N);
pc_m=zeros(1,N);tc_m=zeros(1,N);c=zeros(1,N);mc=zeros(1,N);
for i=1:N
    
    Tr(i)=T1/Tc(i);   
     a(i)=(((0.42747*(R^2)*(Tc(i)^2)))/Pc(i));
     b(i)=(0.08664*R*Tc(i))/Pc(i);
    mb(i)=Yi(i)*b(i);       
    pc_m(i)=Yi(i)*Pc(i);
    tc_m(i)=Yi(i)*Tc(i);
      c1(i)=(-45.7247*((1/3)-Zc(i)));
  c2(i)=((-2.184*exp(c1(i)))+0.2658);
  c(i)=((((1/3)-Zc(i))*((R*Tc(i))/Pc(i)))*c2(i));

   mc(i)=Yi(i)*c(i); 
    
end
% x2=(sum(Tc_m)); %alternative option
% x1=(sum(Pc_m));
% x0=(R*T)/(x1);
b_m=sum(mb);
c_m=sum(mc);
% Tc_m=sum(tc_m);
switch lower(alpha_function)

   case 'original'
        k=zeros(1,N);
        alpha=zeros(1,N);
        first_order_alpha=zeros(1,N);
       for i=1:N
           
        k(i)=0.480+(1.574*w(i))-(0.176*w(i)^2);
        alpha(i)=((1+(k(i)*(1-(Tr(i)^0.5))))^2);
        first_order_alpha(i)=-(k(i)/T1)*(sqrt(Tr(i)))*(sqrt(alpha(i)));
        
       end

        case 'ozokwelu'
        m=zeros(1,N);n2=zeros(1,N);
        alpha=zeros(1,N);
        first_order_alpha=zeros(1,N);
            for i=1:N
            
         m(i)=0.266+(0.4459*w(i)^0.5);
         n2(i)=(1/m(i))*(0.2469+(0.7495*w(i)));

       alpha(i)=((exp(m(i)*(1-((Tr(i))^n2(i)))))^2);
       first_order_alpha(i)=((-2*m(i)*n2(i))/Tc(i))*(Tr(i)^(n2(i)-1))*(exp(2*m(i)*(1-(Tr(i)^n2(i)))));
        
            end

            
                case 'trebble'

           q1=zeros(1,N);
        alpha=zeros(1,N);
        first_order_alpha=zeros(1,N);
        
           for i=1:N
               
         if w(i)<-0.1
             
             q1(i)=0.66208+(4.63961*w(i))+(7.45183*w(i)^2);
             
         elseif w(i)>=-0.1 && w(i)<=0.4
             
             q1(i)=0.35+(0.7924*w(i))+(0.1875*w(i)^2)-(28.93*((0.3-Zc(i))^2));
             
         elseif w(i)>0.4
             
             q1(i)=0.32+(0.9424*w(i))-(28.93*((0.3-Zc(i))^2));
             
         end

         alpha(i)=exp(q1(i)*(1-Tr(i)));
         first_order_alpha(i)=-(q1(i)/Tc(i))*alpha(i);
           
           end
           
           
             case 'nasrifar'
        
            k=zeros(1,N);
            b1=zeros(1,N);b2=zeros(1,N);
            b3=zeros(1,N);first_order_alpha=zeros(1,N);
            alpha=zeros(1,N);second_order_alpha=zeros(1,N);
            
            for i=1:N
                
                if Tr(i)>1
                    
                k(i)=0.480+(1.574*w(i))-(0.175*w(i)^2);
                b1(i)=0.25*(12-(11*k(i))+(k(i)^2));
                b2(i)=0.5*(-6+(9*k(i))-(k(i)^2));
                b3(i)=0.25*(4-(7*k(i))+(k(i)^2));
                
                 alpha(i)=(((b1(i)/(T1/Tc(i)))+(b2(i)/...
        ((T1/Tc(i)))^2)+(b3(i)/((T1/Tc(i)))^3)));
    
    first_order_alpha(i)=-((Tc(i)*b1(i))/(T1^2))-...
        ((2*(Tc(i)^2)*b2(i))/(T1^3))-...
        ((3*(Tc(i)^3)*b3(i))/(T1^4));
     
                elseif Tr(i)<1 || Tr(i)==1                      
            
         
           
        k(i)=0.480+(1.574*w(i))-(0.176*w(i)^2);
        alpha(i)=((1+(k(i)*(1-(Tr(i)^0.5))))^2);
        first_order_alpha(i)=-(k(i)/T1)*(sqrt(Tr(i)))*(sqrt(alpha(i)));
        second_order_alpha(i)=(-k(i)*((-k(i)*sqrt(T1)-...
                            (sqrt(Tc(i))*sqrt(alpha(i))))))/(2*Tc(i)*(T1^1.5));
                        
                end

            end
end

  Cx1=zeros(N,N);
       for i=1:N
    for j=1:N
   
        Cx1(i,j)=((1-Kij(i,j))/2)*(((a(i)*first_order_alpha(i)*...
            a(j)*alpha(j))+(a(j)*first_order_alpha(j)*...
            a(i)*alpha(i)))/(sqrt(((a(i)*alpha(i)*a(j)*alpha(j))))));   
    end
       end
       
F_D_sxV=bsxfun(@times,bsxfun(@times, Yi,Cx1),Yi');
F_D_sxV=sum(F_D_sxV(:));F_D_a_alpha_m=F_D_sxV;


Cx2=zeros(N,N);
for i=1:N
    for j=1:N
   
        Cx2(i,j)=((a(i)*alpha(i)*a(j)*alpha(j))^0.5)*(1-Kij(i,j));
    
    end
end

sxV=bsxfun(@times,bsxfun(@times, Yi,Cx2),Yi');
sxV=sum(sxV(:));a_alpha_m=sxV;


OBJ= -(((((F_D_a_alpha_m)*R*T1)-...
                (R*a_alpha_m)))/...
    ((R*T1)^2))-((b_m-c_m-((a_alpha_m)/(R*T1)))/T1);



end