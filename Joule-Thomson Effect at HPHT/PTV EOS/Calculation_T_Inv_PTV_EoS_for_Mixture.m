function OBJ=Calculation_T_Inv_PTV_EoS_for_Mixture(T,Tc,Pc,w,Kij,Yi,N,Zc)




R=83.14472;
a=zeros(1,N);b=zeros(1,N);Oa=zeros(1,N);mc=zeros(1,N);
alpha=zeros(1,N);first_order_alpha=zeros(1,N);
Tr=zeros(1,N);mb=zeros(1,N);k=zeros(1,N);c=zeros(1,N);
Pc_m=zeros(1,N);Tc_m=zeros(1,N);Ob=zeros(1,N);Oc=zeros(1,N);

for i=1:N
    
    Tr(i)=T/Tc(i);   
          Oa(i)=(0.66121-(0.76105*Zc(i)));
 Ob(i)=(0.02207+(0.20868*Zc(i)));
   Oc(i)=(0.57765-(1.87080*Zc(i)));

    a(i)=(((Oa(i)*(R^2)*(Tc(i)^2)))/Pc(i));
   b(i)=((Ob(i)*R*Tc(i))/Pc(i));
 c(i)=(Oc(i)*R*Tc(i))/Pc(i);
    k(i)=0.46283+(3.58230*w(i)*Zc(i))+(8.19417*(w(i)^2)*(Zc(i)^2));
        alpha(i)=((1+(k(i)*(1-(Tr(i)^0.5))))^2);
        first_order_alpha(i)=-(k(i)/T)*(sqrt(Tr(i)))*(sqrt(alpha(i)));
    mb(i)=Yi(i)*b(i);       
    Pc_m(i)=Yi(i)*Pc(i);
    Tc_m(i)=Yi(i)*Tc(i);
    mc(i)=c(i)*Yi(i);
    
end
% x2=(sum(Tc_m)); %alternative option
% x1=(sum(Pc_m));
% x0=(R*T)/(x1);
b_m=sum(mb);
% c_m=sum(mc);
     
       


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


OBJ= -(((((F_D_a_alpha_m)*R*T)-...
                (R*a_alpha_m)))/...
    ((R*T)^2))-((b_m-((a_alpha_m)/(R*T)))/T);



end