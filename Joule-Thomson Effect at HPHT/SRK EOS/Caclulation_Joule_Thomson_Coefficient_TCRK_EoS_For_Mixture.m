function J_T_m=Caclulation_Joule_Thomson_Coefficient_TCRK_EoS_For_Mixture...
    (Tc,Pc,T,P,Zra_p,M_p,L_p,N_p,N1,Yi,Kij)


% R=83.14472;
  R=8.3144598;
% R=10.732;

a=zeros(1,N1);b=zeros(1,N1);c=zeros(1,N1);
Tr=zeros(1,N1);mb=zeros(1,N1);mc=zeros(1,N1);
Pc_m=zeros(1,N1);Tc_m=zeros(1,N1);
alpha=zeros(1,N1);first_order_alpha=zeros(1,N1);
second_order_alpha=zeros(1,N1);

for i=1:N1
    
    Tr(i)=T/Tc(i);   
   
    a(i)=(0.457235*(R^2)*(Tc(i)^2))/Pc(i);
    c(i)=((R*Tc(i))/Pc(i))*(0.1398-(0.5294*Zra_p(i)));        
    b(i)=((0.07780*R*Tc(i))/Pc(i))-c(i);
    
    mb(i)=Yi(i)*b(i);       
    Pc_m(i)=Yi(i)*Pc(i);
    Tc_m(i)=Yi(i)*Tc(i);
   
    mc(i)=c(i)*Yi(i);
    
     alpha(i)=(Tr(i)^(N_p(i)*(M_p(i)-1)))*(exp(L_p(i)*(1-(Tr(i)^(M_p(i)*N_p(i))))));

    first_order_alpha(i)=(((M_p(i)-1)*(T^((N_p(i)*(M_p(i)-1))-1))*N_p(i)*...
        (exp(L_p(i)*(1-(Tr(i)^(N_p(i)*M_p(i)))))))/(Tc(i)^(N_p(i)*(M_p(i)-1))))-...
    (L_p(i)*N_p(i)*(exp(L_p(i)*(1-(Tr(i)^(N_p(i)*M_p(i))))))*M_p(i)*...
    (T^((2*N_p(i)*M_p(i))-N_p(i)-1)))/(Tc(i)^(N_p(i)*(-1+(2*M_p(i)))));

    second_order_alpha(i)=((N_p(i)*(M_p(i)-1)*((((N_p(i)*(M_p(i)-1))-1)*(T^((N_p(i)*M_p(i))-N_p(i)-2))*...
        (Tc(i)^(N_p(i)*M_p(i)))*(exp(L_p(i)*(1-(Tr(i)^(N_p(i)*M_p(i)))))))-...
    ((T^((2*N_p(i)*M_p(i))-N_p(i)-2))*N_p(i)*L_p(i)*M_p(i)*(exp(L_p(i)*(1-(Tr(i)^(N_p(i)*M_p(i)))))))))/...
    (Tc(i)^(N_p(i)*((2*M_p(i))-1))))-...
    (((M_p(i)*N_p(i)*L_p(i)*((((2*N_p(i)*M_p(i))-N_p(i)-1)*(T^((2*N_p(i)*M_p(i))-N_p(i)-2))*(Tc(i)^(M_p(i)*N_p(i)))*...
    (exp(L_p(i)*(1-(Tr(i)^(N_p(i)*M_p(i)))))))-...
    (M_p(i)*N_p(i)*L_p(i)*(T^((3*N_p(i)*M_p(i))-N_p(i)-2))*(exp(L_p(i)*(1-(Tr(i)^(N_p(i)*M_p(i)))))))))/(Tc(i)^(N_p(i)*((3*M_p(i))-1)))));
    
    
end
% x2=(sum(Tc_m)); %alternative option
% x1=(sum(Pc_m));
% x0=(R*T)/(x1);
b_m=sum(mb);
c_m=sum(mc);
% Tc_m=sum(tc_m);
    x0=((R*(T))/P);


         Cx1=zeros(N1,N1);Cx2=zeros(N1,N1);Cx3=zeros(N1,N1);

       for i=1:N1
    for j=1:N1
   
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
        
                
        F=@(x) ((R*T)/(x-b_m))-((a_alpha_m)/(((x+c_m)*(x+b_m+(2*c_m)))+((b_m+c_m)*(x-b_m))))-P;
         
        V=fzero(F,x0);
        
        
        
        cp_star_m=zeros(1,N1);cp_star=zeros(1,N1);
        for i=1:N1
        
         cp_star(i)=A(i)+(B(i)*T)+(C(i)*(T^2))+(D(i)*(T^3))+(E(i)*(T^4));
         cp_star_m(i)=cp_star(i)*Yi(i);
         
        end
         Cv_star_m=sum(cp_star_m)-R;
    
%     Tr=T/Tc;
    
%     Z=(P*OF)/(R*T);
                   

          Cv_m=Cv_star_m-(((T*S_D_a_alpha_m)/(2*sqrt(2)*(b_m+c_m)))*...
          (log((V-(sqrt(2)*(b_m+c_m))+(2*c_m)+b_m)/(V+(sqrt(2)*(b_m+c_m))+(2*c_m)+b_m))));
     
     
     first_order_P_T_m=(R/(V-b_m))-((F_D_a_alpha_m)/(((V+c_m)*(V+b_m+(2*c_m)))+((b_m+c_m)*(V-b_m))));
    
     first_order_P_V_m=((-R*T)/((V-b_m)^2))+((2*a_alpha_m*(V+(2*c_m)+b_m))/(((V^2)+(4*c_m*V)+(2*b_m*V)+(2*(c_m^2))-(b_m^2))^2));
     
    
    Cp_m=Cv_m-((T*(first_order_P_T_m^2))/(first_order_P_V_m));
    
     
    J_T_m=((T*((-first_order_P_T_m)/first_order_P_V_m))-V)/Cp_m;
        
        
        
        
        
end
    
