clc
clear 
close all

format long

        data=load('SS1');
        
        Tc=data.Tc;
Pc=data.Pc;
w=data.w;

n1=numel(Tc);
P=1e-7;R=83.14472;
k=zeros(n1,1);
b=zeros(n1,1);
a=zeros(n1,1);
for i=1:n1
    
  a(i)=(0.45724*(R^2)*(Tc(i)^2))/Pc(i);
  b(i)=(0.07780*R*Tc(i))/Pc(i);
  k(i)=0.37464+(1.54226*w(i))-(0.26992*w(i)^2);

end


x1=input('Enter lowe range of temperature, T1=');
x2=input('Enter uper range of temperature, T2=');


 T=linspace(x1,x2,2500)';n=numel(T);
 
 
 x0=zeros(n,1);
for i=1:n
    
    
    x0(i)=(R*T(i))/P;
    
    
end


OF=zeros(n,n1);
for j=1:n1
    
for i=1:n
    
     f=@(x) ((R*T(i))/(x-b(j)))-((a(j)*((1+(k(j)*(1-((T(i)/Tc(j))^0.5))))^2))/((x^2)+(2*b(j)*x)-(b(j)^2)))-P;
     
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
    Z(i,j)=(P*OF(i,j))/(R*T(i));
    
    alpha(i,j)=((1+(k(j)*(1-((T(i)/Tc(j))^0.5))))^2); % Original
    A(i,j)=(a(j)*alpha(i,j))/((R*T(i))^2);
    
    
    
    DEI(i,j)=((a(j)*(alpha(i,j)+(k(j)*sqrt(Tr(i,j))*...
        sqrt(alpha(i,j)))))); % Original
    
    
    Delta_Enthalpy(i,j)=((((DEI(i,j))/(2*sqrt(2)*b(j))))*...
        (log((OF(i,j)+((1-sqrt(2))*b(j)))/(OF(i,j)+((1+sqrt(2))*b(j))))))+(R*T(i)*(Z(i,j)-1)); % Original
    
    Delta_Entropy(i,j)=(R*(log(Z(i,j)-B(i,j))))+...
        ((((a(j)*k(j)*sqrt(Tr(i,j))*sqrt(alpha(i,j)))/T(i))/(2*sqrt(2)*b(j)))*...
        (log((OF(i,j)+((1-sqrt(2))*b(j)))/(OF(i,j)+((1+sqrt(2))*b(j)))))); % Original
    
    
    Delta_Helmholtz_Energy(i,j)=(R*T(i)*(log(1/(Z(i,j)-B(i,j)))))-...
        (((a(j)*alpha(i,j))/(2*sqrt(2)*b(j)))*...
        (log((OF(i,j)+((1-sqrt(2))*b(j)))/(OF(i,j)+((1+sqrt(2))*b(j))))));
    
    
    Delta_Internal_Energy(i,j)=((((DEI(i,j))/(2*sqrt(2)*b(j))))*...
        (log((OF(i,j)+((1-sqrt(2))*b(j)))/(OF(i,j)+((1+sqrt(2))*b(j))))));

    
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

% I=max(F_B);
% disp('maximum of Residual Enthalpy per Residual Entropy=');disp(I);
% 
% Andis=find(F_B==I);
% % [row,colm]=find(OBV==I)
% 
% T_Final=T(Andis);
% disp('The Boyle temperature, T');disp(T_Final)
% 
% RH=Delta_Enthalpy(Andis);
% RS=Delta_Entropy(Andis);
% 
% 
% plot(T_Final,I,'ok',T,F_B)
% 
% % legend({'Boyle Temperature','Ratio of Residual Enthalpy to Residual Entropy'})
