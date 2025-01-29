clc
clear 
close all

format long

        data=load('SS1');
        
Tc=data.Tc;Pc=data.Pc;
w=data.w;Zc=data.Zc;

n1=numel(Tc);
P=1e-7;R=83.14472;
F=zeros(n1,1);b=zeros(n1,1);a=zeros(n1,1);c=zeros(n1,1);
Oa=zeros(n1,1);Ob=zeros(n1,1);Oc=zeros(n1,1);
for i=1:n1
    

   Oa(i)=(0.66121-(0.76105*Zc(i)));
   Ob(i)=(0.02207+(0.20868*Zc(i)));
   Oc(i)=(0.57765-(1.87080*Zc(i)));

    a(i)=(((Oa(i)*(R^2)*(Tc(i)^2)))/Pc(i));
   b(i)=((Ob(i)*R*Tc(i))/Pc(i));
 c(i)=(Oc(i)*R*Tc(i))/Pc(i);
    F(i)=0.46283+(3.58230*w(i)*Zc(i))+(8.19417*(w(i)^2)*(Zc(i)^2));
    
end


x1=input('Enter lowe range of temperature, T1=');
x2=input('Enter uper range of temperature, T2=');


 T=linspace(x1,x2,5000)';n=numel(T);
 
 
 x0=zeros(n,1);
for i=1:n
    
    
    x0(i)=(R*T(i))/P;
    
    
end


OF=zeros(n,n1);
for j=1:n1
    
for i=1:n
    
     f=@(x) ((R*T(i))/(x-b(j)))-((a(j)*((1+(F(j)*(1-((T(i)/Tc(j))^0.5))))^2))/((x*(x+b(j)))+(c(j)*(x-b(j)))))-P;
     
     OF(i,j)=fzero(f,x0(i));
%           OF(i)=fsolve(f,x0(i));
% OF(i)=fminsearch(f,x0(i))

end

end

DEI=zeros(n,1);C=zeros(n,1);
A=zeros(n,1);alpha=zeros(n,1);B=zeros(n,1);Tr=zeros(n,1);
Delta_Gipss=zeros(n,1);Delta_Helmholtz_Energy=zeros(n,1);
Delta_Internal_Energy=zeros(n,1);Delta_Entropy=zeros(n,1);
Delta_Enthalpy=zeros(n,1);Z=zeros(n,1);F_B=zeros(n,1);



for j=1:n1

for i=1:n
    
    Tr(i,j)=T(i)/Tc(j);
    B(i,j)=(b(j)*P)/(R*T(i));
    C(i,j)=(c(j)*P)/(R*T(i));
    Z(i,j)=(P*OF(i,j))/(R*T(i));
    alpha(i,j)=((1+(F(j)*(1-((T(i)/Tc(j))^0.5))))^2);
    
    A(i,j)=(a(j)*alpha(i,j))/((R*T(i))^2);
     
    DEI(i,j)=(a(j)*alpha(i,j)-(T(i)*(((-F(j)*a(j))/T(i))*sqrt(Tr(i,j))*sqrt(alpha(i,j)))));
  
    Delta_Enthalpy(i,j)=(R*T(i)*(Z(i,j)-1))+...
        ((DEI(i,j)/(sqrt((c(j)^2)+(6*b(j)*c(j))+(b(j)^2))))*(log((((2*OF(i,j))-...
        (sqrt((c(j)^2)+(6*b(j)*c(j))+(b(j)^2)))+c(j)+b(j))/((2*OF(i,j))+...
        (sqrt((c(j)^2)+(6*b(j)*c(j))+(b(j)^2)))+c(j)+b(j))))));
    
    Delta_Entropy(i,j)=(R*(log(Z(i,j)-B(i,j))))-...
        (((((-F(j)*a(j))/T(i))*sqrt(Tr(i,j))*sqrt(alpha(i,j))/...
        (sqrt((c(j)^2)+(6*b(j)*c(j))+(b(j)^2)))))*(log((((2*OF(i,j))-...
        (sqrt((c(j)^2)+(6*b(j)*c(j))+(b(j)^2)))+c(j)+b(j))/((2*OF(i,j))+...
        (sqrt((c(j)^2)+(6*b(j)*c(j))+(b(j)^2)))+c(j)+b(j))))));
    
        Delta_Internal_Energy(i,j)=Delta_Enthalpy(i,j)-(R*T(i)*(Z(i,j)-1));
    
    Delta_Helmholtz_Energy(i,j)=((R*T(i))*(log(1/(Z(i,j)-B(i,j)))))-...
        (((a(j)*alpha(i,j))/(sqrt((c(j)^2)+(6*b(j)*c(j))+(b(j)^2))))*(log((((2*OF(i,j))-...
        (sqrt((c(j)^2)+(6*b(j)*c(j))+(b(j)^2)))+c(j)+b(j))/((2*OF(i,j))+...
        (sqrt((c(j)^2)+(6*b(j)*c(j))+(b(j)^2)))+c(j)+b(j))))));

    F_B(i,j)=Delta_Enthalpy(i,j)/Delta_Entropy(i,j);
    
end


end

I=zeros(n1,1);
for i=1:n1
    
  I(i)=max(F_B(:,i));

end

Andis=zeros(n1,1);
for i=1:n1
    
  Andis(i)=find(F_B(:,i)==I(i));

end

T_Final=zeros(n1,1);
for i=1:n1

    T_Final(i)=T(Andis(i));
    
end

