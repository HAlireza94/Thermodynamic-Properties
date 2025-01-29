clc
clear
close all
format long

data=load('SS3');

 Tc=data.Tc;Pc=data.Pc;
        w=data.w;Zra=data.Zra;
        L=data.L;M=data.M;N=data.N;
        n1=numel(Tc);

R=83.14472;

c=zeros(n1,1);
b=zeros(n1,1);
a=zeros(n1,1);
for i=1:n1
    
     a(i)=(1/(9*((2^(1/3))-1)))*((((R^2)*(Tc(i)^2)))/Pc(i));
            c(i)=((R*Tc(i))/Pc(i))*(0.2150-(0.7314*Zra(i)));
            b(i)=((((2^(1/3))-1)/3)*((R*Tc(i))/Pc(i)))-c(i);

end

T=linspace(200,2500,2000000)';nt=numel(T);

ax1=zeros(nt,n1);
B1=zeros(nt,n1);

Tr=zeros(nt,n1);alpha=zeros(nt,n1);
for j=1:n1
for i=1:nt
    % Orginal 1976
    Tr(i,j)=T(i)/Tc(j);
                   alpha(i,j)=(Tr(i,j)^(N(j)*(M(j)-1)))*(exp(L(j)*(1-(Tr(i,j)^(M(j)*N(j))))));

    B1(i,j)=b(j)-((a(j)*alpha(i,j))/(R*T(i)));
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