function P=Calculation_JTIC_VTSRK_EoS_For_Mixture(alpha_function,T,Tc,Pc,w,N,Yi,Kij,Zc)


R=83.14472;
a=zeros(1,N);b=zeros(1,N);c1=zeros(1,N);
Tr=zeros(1,N);mb=zeros(1,N);c2=zeros(1,N);
Pc_m=zeros(1,N);Tc_m=zeros(1,N);c=zeros(1,N);mc=zeros(1,N);
for i=1:N
    
    Tr(i)=T/Tc(i);   
     a(i)=(((0.42747*(R^2)*(Tc(i)^2)))/Pc(i));
     b(i)=(0.08664*R*Tc(i))/Pc(i);
    mb(i)=Yi(i)*b(i);       
    Pc_m(i)=Yi(i)*Pc(i);
    Tc_m(i)=Yi(i)*Tc(i);
    c1(i)=(-45.7247*((1/3)-Zc(i)));
  c2(i)=((-2.184*exp(c1(i)))+0.2658);
  c(i)=((((1/3)-Zc(i))*((R*Tc(i))/Pc(i)))*c2(i));

   mc(i)=Yi(i)*c(i);        
    
end
x2=(sum(Tc_m)); %alternative option
x1=(sum(Pc_m));
x0=(R*(T-(0.01*x2)))/(x1);
b_m=sum(mb);
c_m=sum(mc);
switch lower(alpha_function)

   case 'original'
        k=zeros(1,N);
        alpha=zeros(1,N);
        first_order_alpha=zeros(1,N);
       for i=1:N
           
        k(i)=0.480+(1.574*w(i))-(0.176*w(i)^2);
        alpha(i)=((1+(k(i)*(1-(Tr(i)^0.5))))^2);
        first_order_alpha(i)=-(k(i)/T)*(sqrt(Tr(i)))*(sqrt(alpha(i)));
        
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
                
                 alpha(i)=(((b1(i)/(T/Tc(i)))+(b2(i)/...
        ((T/Tc(i)))^2)+(b3(i)/((T/Tc(i)))^3)));
    
    first_order_alpha(i)=-((Tc(i)*b1(i))/(T^2))-...
        ((2*(Tc(i)^2)*b2(i))/(T^3))-...
        ((3*(Tc(i)^3)*b3(i))/(T^4));
     

    second_order_alpha(i)=((2*Tc(i)*b1(i))/(T^3))+...
        ((6*(Tc(i)^2)*b2(i))/(T^4))+...
        ((12*(Tc(i)^3)*b3(i))/(T^5));

                elseif Tr(i)<1 || Tr(i)==1                      
            
         
           
        k(i)=0.480+(1.574*w(i))-(0.176*w(i)^2);
        alpha(i)=((1+(k(i)*(1-(Tr(i)^0.5))))^2);
        first_order_alpha(i)=-(k(i)/T)*(sqrt(Tr(i)))*(sqrt(alpha(i)));
        second_order_alpha(i)=(-k(i)*((-k(i)*sqrt(T)-...
                            (sqrt(Tc(i))*sqrt(alpha(i))))))/(2*Tc(i)*(T^1.5));
                        
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

             
    
    F1=@(x) (T*((R/(x+c_m-b_m))-(F_D_a_alpha_m/((x+c_m)*(x+c_m+b_m)))))+...
        (x*(((-R*T)/((x+c_m-b_m)^2))+...
        ((a_alpha_m*((2*c_m)+(2*x)+b_m)))/((((x+c_m)*(x+c_m+b_m)))^2)));
    
    V=fzero(F1,x0);
    
    P=((R*T)/(V+c_m-b_m))-((a_alpha_m)/((V+c_m)*(V+c_m+b_m))); %$this was wrong
                
end
    
