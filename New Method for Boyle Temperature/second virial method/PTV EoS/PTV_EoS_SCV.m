clc
clear
close all
format long

data=load('SS1');

Tc=data.Tc;
Pc=data.Pc;
w=data.w;
Zc=data.Zc;
n1=numel(Tc);

R=83.14472;

k=zeros(n1,1);b=zeros(n1,1);a=zeros(n1,1);c=zeros(n1,1);
Oa=zeros(n1,1);Ob=zeros(n1,1);Oc=zeros(n1,1);
for i=1:n1
    
  Oa(i)=(0.66121-(0.76105*Zc(i)));
   Ob(i)=(0.02207+(0.20868*Zc(i)));
   Oc(i)=(0.57765-(1.87080*Zc(i)));

    a(i)=(((Oa(i)*(R^2)*(Tc(i)^2)))/Pc(i));
   b(i)=((Ob(i)*R*Tc(i))/Pc(i));
 c(i)=(Oc(i)*R*Tc(i))/Pc(i);
    k(i)=0.46283+(3.58230*w(i)*Zc(i))+(8.19417*(w(i)^2)*(Zc(i)^2));

end

T=linspace(200,2500,5000000)';nt=numel(T);

ax1=zeros(nt,n1);
B1=zeros(nt,n1);

Tr=zeros(nt,n1);
for j=1:n1
for i=1:nt
    % Orginal 1976
    Tr(i,j)=T(i)/Tc(j);
    
    B1(i,j)=b(j)-((a(j)*((1+(k(j)*(1-(T(i)/Tc(j))^0.5)))^2))/(R*T(i)));
    ax1(i,j)=((B1(i)*Pc(j))/(R*Tc(j)));
     
end
end

I=zeros(n1,1);
for i=1:n1
    
  I(i)=min(abs(B1(:,i)));

end

Andis=zeros(n1,1);
for i=1:n1
    
  Andis(i)=find(abs(B1(:,i))==I(i));

end

T_Final=zeros(n1,1);
for i=1:n1

    T_Final(i)=T(Andis(i));
    
end

% T_Final=zeros(n1,1);
% for i=1:n1
%     
%     T_Final(i)=Tr_Final(i)*Tc(i);
%     
% end