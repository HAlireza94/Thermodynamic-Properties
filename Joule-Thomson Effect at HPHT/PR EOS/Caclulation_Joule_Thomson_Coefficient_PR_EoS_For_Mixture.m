function [J_T_m, V]=Caclulation_Joule_Thomson_Coefficient_PR_EoS_For_Mixture...
    (alpha_function,T,Tc,Pc,w,N,Yi,Kij,A,B,C,D,E,P)


% R=83.14472;
  R=8.3144598;

a=zeros(1,N);b=zeros(1,N);
Tr=zeros(1,N);mb=zeros(1,N);
Pc_m=zeros(1,N);tc_m=zeros(1,N);
for i=1:N
    
    Tr(i)=T/Tc(i);   
    a(i)=(0.457235*(R^2)*(Tc(i)^2))/Pc(i);
    b(i)=(0.077796*R*Tc(i))/Pc(i);
    mb(i)=Yi(i)*b(i);       
    Pc_m(i)=Yi(i)*Pc(i);
    tc_m(i)=Yi(i)*Tc(i);
    
    
end
% x2=(sum(Tc_m)); %alternative option
x1=(sum(Pc_m));
% x0=(R*T)/(x1);
b_m=sum(mb);
Tc_m=sum(tc_m);
    x0=((R*(T+Tc_m+200))/(P));

switch lower(alpha_function)

   case 'original'
        k=zeros(1,N);
        alpha=zeros(1,N);second_order_alpha=zeros(1,N);
        first_order_alpha=zeros(1,N);
       for i=1:N
           
        k(i)=0.37464+(1.54226*w(i))-(0.26992*w(i)^2);
        alpha(i)=((1+(k(i)*(1-(Tr(i)^0.5))))^2);
        first_order_alpha(i)=-(k(i)/T)*(sqrt(Tr(i)))*(sqrt(alpha(i)));
        second_order_alpha(i)=(-k(i)*((-k(i)*sqrt(T)-...
                            (sqrt(Tc(i))*sqrt(alpha(i))))))/(2*Tc(i)*(T^1.5));
        
       end

        case 'coquelet'
%         k1=zeros(1,N);
%         alpha=zeros(1,N);second_order_alpha=zeros(1,N);
%         first_order_alpha=zeros(1,N);
            for i=1:N
      if Tr(i)>1
          
        k1(i)=0.41287+(1.34494*w(i))+(0.00421*w(i)^2);
        alpha(i)=exp(k1(i)*(1-Tr(i)));
        first_order_alpha(i)=-(k1(i)/Tc(i))*alpha(i);
        second_order_alpha(i)=((k1(i)/Tc(i))^2)*alpha(i);
        
      elseif Tr(i)<1
          k1(i)=0.41287+(1.34494*w(i))+(0.00421*w(i)^2);
       k2(i)=0.03982+(0.08551*w(i))-(1.05521*w(i)^2);
            k3(i)=-0.01871+(1.93328*w(i));
       alpha(i)=(((exp(k1(i)*(1-Tr(i)))))*((1+(k2(i)*...
           (1-sqrt(Tr(i)))^2)+(k3(i)*(1-sqrt(Tr(i)))^3))^2));
     
        x1(i)=(1/(Tc(i)^4))*((((Tc(i)^1.5)*(1+k2(i)+k3(i)-(k3(i)*Tr(i)^1.5)))-...
    (Tc(i)*sqrt(T)*((2*k2(i))+(3*k3(i))))+(T*sqrt(Tc(i))*(k2(i)+(3*k3(i)))))^2);

x2(i)=1+k2(i)+k3(i)-(sqrt(Tr(i))*((2*k2(i))+...
    (3*k3(i))))+(Tr(i)*(k2(i)+(3*k3(i))))-(k3(i)*Tr(i)^1.5);

x3(i)=((-((2*Tc(i)*k2(i)*(1-sqrt(Tr(i))))+((3*k3(i)*...
    (sqrt(Tc(i))-sqrt(T)))^2)))/(2*sqrt(T)*Tc(i)^1.5));

y1(i)=((k2(i)+(3*k3(i)))/Tc(i))-(((2*k2(i))+(3*k3(i)))/...
    (2*sqrt(Tc(i))*sqrt(T)))-((1.5*k3(i)*sqrt(Tr(i)))/Tc(i));

y2(i)=x3(i);

y3(i)=(-((2*sqrt(Tc(i))*k2(i)*(-sqrt(T)-(sqrt(Tc(i))*...
    (1-sqrt(Tr(i))))))+(3*k3(i)*(T-Tc(i))))/(4*(Tc(i)*T)^1.5));

y4(i)=x2(i);

x4(i)=((k1(i)/Tc(i))^2)*(((exp(k1(i)*(1-Tr(i)))))*...
    ((1+(k2(i)*(1-sqrt(Tr(i)))^2)+(k3(i)*(1-sqrt(Tr(i)))^3))^2));

x5(i)=2*(exp(k1(i)*(1-Tr(i))))*((y1(i)*y2(i))+(y3(i)*y4(i)));

x6(i)=((-4*k1(i))/Tc(i))*y2(i)*y4(i)*(exp(k1(i)*(1-Tr(i))));

first_order_alpha(i)=(exp(k1(i)*(1-Tr(i))))*((2*x2(i)*x3(i))-(k1(i)*x1(i)));
second_order_alpha(i)=x4(i)+x5(i)+x6(i);
          
      end

            end
            
        case 'haghtalab'
            
        n2=zeros(1,N);
        alpha=zeros(1,N);second_order_alpha=zeros(1,N);
        first_order_alpha=zeros(1,N);
            for i=1:N
             
            n2(i)=1.7309+(1.6571*w(i))+(0.1554*w(i)^2);
            alpha(i)=(exp(1-(n2(i)^(log(Tr(i))))));    
            first_order_alpha(i)=((-(exp(1-(n2(i)^(log(Tr(i))))))*...
                (T^(log(n2(i))-1))*(log(n2(i))))/(Tc(i)^(log(n2(i)))));

            second_order_alpha(i)=(((-log(n2(i)))*(exp(1-(n2(i)^(log(Tr(i))))))*...
            (Tr(i)^(log(n2(i))))*((log(n2(i)))-1-((log(n2(i)))*(Tr(i)^(log(n2(i)))))))/(T^2));
            
            end
            
        case 'saffari'
           
            k1=zeros(1,N);
            k2=zeros(1,N);y1=zeros(1,N);
            k3=zeros(1,N);y2=zeros(1,N);
            alpha=zeros(1,N);second_order_alpha=zeros(1,N);
            first_order_alpha=zeros(1,N);
            
            for i=1:N
                
            k1(i)=0.003091+(0.013145*w(i));
            k2(i)=-0.006478+(0.482173*w(i));
            k3(i)=(3.58616*w(i))+0.721306;
            alpha(i)=(exp((k1(i)*Tr(i))+(k2(i)*(log(Tr(i))))+(k3(i)*(1-(Tr(i)^0.5)))));
            
            first_order_alpha(i)=(exp((k1(i)*Tr(i))+(k3(i)*(1-(Tr(i)^0.5)))))*...
                (((2*k2(i)*(Tc(i)^1.5)*(T^(k2(i)-1)))-(k3(i)*Tc(i)*(T^(k2(i)-0.5)))+...
                (2*k1(i)*sqrt(Tc(i))*(T^k2(i))))/(2*(Tc(i)^(k2(i)+1.5))));
            
            y1(i)=((k2(i)*(T^(k2(i)-1)))/(Tc(i)^k2(i)))*((k2(i)/T)-...
                (k3(i)/(2*sqrt(Tc(i))*sqrt(T)))+(k1(i)/Tc(i)));
            
            y2(i)=(Tr(i)^k2(i))*((k3(i)/(4*(T^1.5)*sqrt(Tc(i))))-(k2(i)/(T^2)));
                                    
            second_order_alpha(i)=((exp((k1(i)*Tr(i))+(k3(i)*(1-(Tr(i)^0.5)))))*...
                ((k1(i)/Tc(i))-(k3(i)/(2*sqrt(Tc(i))*sqrt(T))))*(Tr(i)^k2(i))*...
                ((k2(i)/T)-(k3(i)/(2*sqrt(Tc(i))*sqrt(T)))+(k1(i)/Tc(i))))+...
                ((y1(i)+y2(i))*(exp((k1(i)*Tr(i))+(k3(i)*(1-(Tr(i)^0.5))))));                                                                                                    

            
            end
            
end



         Cx1=zeros(N,N);Cx2=zeros(N,N);Cx3=zeros(N,N);

       for i=1:N
    for j=1:N
   
        Cx1(i,j)=((1-Kij(i,j))/2)*(((a(i)*first_order_alpha(i)*...
            a(j)*alpha(j))+(a(j)*first_order_alpha(j)*...
            a(i)*alpha(i)))/(sqrt(((a(i)*alpha(i)*a(j)*alpha(j))))));   
        
         Cx2(i,j)=((a(i)*alpha(i)*a(j)*alpha(j))^0.5)*(1-Kij(i,j));
         
         Cx3(i,j)=((1-Kij(i,j))/4)*(((((2*a(i)*alpha(i)*a(j)*alpha(j)*((a(i)*second_order_alpha(i)*a(j)*alpha(j))+...
             (2*a(i)*first_order_alpha(i)*a(j)*first_order_alpha(j))+...
             (a(i)*alpha(i)*a(j)*second_order_alpha(j)))))-...
         (((a(i)*first_order_alpha(i)*a(j)*alpha(j))+(a(j)*first_order_alpha(j)*a(i)*alpha(i)))^2))/...
         ((a(i)*alpha(i)*a(j)*alpha(j))^1.5)));     
         
    end
       end
       
F_D_sxV=bsxfun(@times,bsxfun(@times, Yi,Cx1),Yi');
F_D_sxV=sum(F_D_sxV(:));F_D_a_alpha_m=F_D_sxV;

S_D_sxV=bsxfun(@times,bsxfun(@times, Yi,Cx3),Yi');
S_D_sxV=sum(S_D_sxV(:));S_D_a_alpha_m=S_D_sxV;


sxV=bsxfun(@times,bsxfun(@times, Yi,Cx2),Yi');
sxV=sum(sxV(:));a_alpha_m=sxV;

    
    F=@(x) ((R*T)/(x-b_m))-((a_alpha_m)/...
             ((x^2)+(2*b_m*x)-(b_m^2)))-P;
                
        V=fzero(F,x0);
        
        
        
        cp_star_m=zeros(1,N);cp_star=zeros(1,N);
        for i=1:N
        
         cp_star(i)=A(i)+(B(i)*T)+(C(i)*(T^2))+(D(i)*(T^3))+(E(i)*(T^4));
         cp_star_m(i)=cp_star(i)*Yi(i);
         
        end
         Cv_star_m=sum(cp_star_m)-R;
    
%     Tr=T/Tc;
    
%     Z=(P*OF)/(R*T);
              
     Cv_m=Cv_star_m+(T*...
     ((S_D_a_alpha_m)/(2*sqrt(2)*b_m))*(log((V+...
     ((1+sqrt(2))*b_m))/(V+((1-sqrt(2))*b_m)))));
 
    first_order_P_T_m=(R/(V-b_m))-((F_D_a_alpha_m)/...
        ((V^2)+(2*b_m*V)-(b_m^2)));
    
    first_order_P_V_m=(((2*a_alpha_m*(V+b_m))...
        /(((V^2)+(2*b_m*V)-(b_m^2))^2)))-...
        ((R*T)/((V-b_m)^2));
    
    Cp_m=Cv_m-((T*(first_order_P_T_m^2))/(first_order_P_V_m));
    
     
    J_T_m=((T*((-first_order_P_T_m)/first_order_P_V_m))-V)/Cp_m;
        
        
        
        
        
end
    
