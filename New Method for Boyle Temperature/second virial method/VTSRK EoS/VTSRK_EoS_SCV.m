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

k=zeros(n1,1);
b=zeros(n1,1);
a=zeros(n1,1);
c=zeros(n1,1);c1=zeros(n1,1);c2=zeros(n1,1);
for i=1:n1
    
 a(i)=(((0.42747*(R^2)*(Tc(i)^2)))/Pc(i));
 b(i)=(0.08664*R*Tc(i))/Pc(i);

 c1(i)=(-45.7247*((1/3)-Zc(i)));
 c2(i)=((-2.184*exp(c1(i)))+0.2658);
 c(i)=((((1/3)-Zc(i))*((R*Tc(i))/Pc(i)))*c2(i));
 
 k(i)=0.480+(1.574*w(i))-(0.175*w(i)^2);

end

T=linspace(200,2500,5000000)';nt=numel(T);

ax1=zeros(nt,n1);
B1=zeros(nt,n1);

Tr=zeros(nt,n1);
for j=1:n1
for i=1:nt
    % Orginal 1976
    Tr(i,j)=T(i)/Tc(j);
    
    B1(i,j)=b(j)-c(j)-((a(j)*((1+(k(j)*(1-(T(i)/Tc(j))^0.5)))^2))/(R*T(i)));
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