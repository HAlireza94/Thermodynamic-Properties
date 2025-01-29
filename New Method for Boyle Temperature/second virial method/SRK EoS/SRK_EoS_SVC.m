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
for i=1:n1
    
 a(i)=(((0.42747*(R^2)*(Tc(i)^2)))/Pc(i));
 b(i)=(0.08664*R*Tc(i))/Pc(i);
 k(i)=0.480+(1.574*w(i))-(0.175*w(i)^2);

end

T=linspace(200,2500,500000)';nt=numel(T);

ax1=zeros(nt,n1);
B1=zeros(nt,n1);

Tr=zeros(nt,n1);
for j=1:n1
for i=1:nt
    % Orginal 1972
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