function [J_T_m, V]=Caclulation_Joule_Thomson_Coefficient_SRK_EoS_For_Mixture...
    (alpha_function,T,Tc,Pc,w,N,Yi,Kij,A,B,C,D,E,P,Zc)


% R=83.14472;
  R=8.3144598;

a=zeros(1,N);b=zeros(1,N);
Tr=zeros(1,N);mb=zeros(1,N);
pc_m=zeros(1,N);Tc_m=zeros(1,N);
for i=1:N
    
    Tr(i)=T/Tc(i);   
     a(i)=(((0.42747*(R^2)*(Tc(i)^2)))/Pc(i));
     b(i)=(0.08664*R*Tc(i))/Pc(i);
    mb(i)=Yi(i)*b(i);       
    pc_m(i)=Yi(i)*Pc(i);
    Tc_m(i)=Yi(i)*Tc(i);
    
    
    
end
x2=(sum(Tc_m)); %alternative option
% x1=(sum(Pc_m));
% x0=(R*T)/(x1);

b_m=sum(mb);
% Tc_m=sum(tc_m);
x0=((R*(T+500))/P);
% Tr_m=T/Tc_m;+Tc_m+100

switch lower(alpha_function)

   case 'original'
        k=zeros(1,N);
        alpha=zeros(1,N);second_order_alpha=zeros(1,N);
        first_order_alpha=zeros(1,N);
       for i=1:N
           
        k(i)=0.480+(1.574*w(i))-(0.176*w(i)^2);
        alpha(i)=((1+(k(i)*(1-(Tr(i)^0.5))))^2);
        first_order_alpha(i)=-(k(i)/T)*(sqrt(Tr(i)))*(sqrt(alpha(i)));
        second_order_alpha(i)=(-k(i)*((-k(i)*sqrt(T)-...
                            (sqrt(Tc(i))*sqrt(alpha(i))))))/(2*Tc(i)*(T^1.5));

        
       end
       
       case 'jub'
        k=zeros(1,N);
        alpha=zeros(1,N);second_order_alpha=zeros(1,N);
        first_order_alpha=zeros(1,N);
       for i=1:N
           
         k(i)=0.4810+(1.5963*w(i))-(0.2963*w(i)^2)+(0.1223*w(i)^3);
        alpha(i)=((1+(k(i)*(1-(Tr(i)^0.5))))^2);
        first_order_alpha(i)=-(k(i)/T)*(sqrt(Tr(i)))*(sqrt(alpha(i)));
        second_order_alpha(i)=(-k(i)*((-k(i)*sqrt(T)-...
                            (sqrt(Tc(i))*sqrt(alpha(i))))))/(2*Tc(i)*(T^1.5));

        
       end


        case 'ozokwelu'
        m=zeros(1,N);n2=zeros(1,N);
        alpha=zeros(1,N);second_order_alpha=zeros(1,N);
        first_order_alpha=zeros(1,N);
            for i=1:N
            
         m(i)=0.266+(0.4459*w(i)^0.5);
         n2(i)=(1/m(i))*(0.2469+(0.7495*w(i)));

       alpha(i)=((exp(m(i)*(1-((Tr(i))^n2(i)))))^2);
       first_order_alpha(i)=((-2*m(i)*n2(i))/Tc(i))*(Tr(i)^(n2(i)-1))*(exp(2*m(i)*(1-(Tr(i)^n2(i)))));
         second_order_alpha(i)=((-2*m(i)*n2(i)*(T^(n2(i)-2))*((((Tc(i)^n2(i))*(n2(i)-1))-...
                           (2*m(i)*n2(i)*T^n2(i))))*(exp(2*m(i)*(1-(Tr(i)^n2(i))))))/(Tc(i)^(2*n2(i))));
            end
             case 'trebble'

           q1=zeros(1,N);
        alpha=zeros(1,N);second_order_alpha=zeros(1,N);
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
         second_order_alpha(i)=((q1(i)/Tc(i))^2)*alpha(i);

           
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
             (x*(x+b_m)))-P;
                
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
     ((S_D_a_alpha_m)/(b_m))*...
     (log((b_m/V)+1)));
    
     first_order_P_T_m=(R/(V-b_m))-...
         ((F_D_a_alpha_m)/(V*(V+b_m)));
    
     first_order_P_V_m=((-R*T)/((V-b_m)^2))+...
         (a_alpha_m*((2*V)+b_m))/((V^2)*((V+b_m)^2));
    
    Cp_m=Cv_m-((T*(first_order_P_T_m^2))/(first_order_P_V_m));
    
     
    J_T_m=((T*((-first_order_P_T_m)/first_order_P_V_m))-V)/Cp_m;
        
        
        
        
        
end
    
