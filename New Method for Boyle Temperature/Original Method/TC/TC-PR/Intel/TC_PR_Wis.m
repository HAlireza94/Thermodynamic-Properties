clc
clear
close all
format long
P=1e-07;R=83.14472;

data=load('SS2');

 Tc=data.Tc;Pc=data.Pc;
        w=data.w;Zra=data.Zra;
        L=data.L;M=data.M;N=data.N;

n1=numel(Tc);

c=zeros(n1,1);
b=zeros(n1,1);
a=zeros(n1,1);
for i=1:n1
    
   a(i)=(0.45724*(R^2)*(Tc(i)^2))/Pc(i);
            c(i)=((R*Tc(i))/Pc(i))*(0.1975-(0.7325*Zra(i)));
            b(i)=((0.07780*R*Tc(i))/Pc(i))-c(i);

end


x1=input('Enter lowe range of temperature, T1=');
x2=input('Enter uper range of temperature, T2=');


 T=linspace(x1,x2,1000)';n=numel(T);

 x0=zeros(n,1);
for i=1:n
    
    
    x0(i)=(R*T(i))/P;
    
    
end


OF=zeros(n,n1);Tr=zeros(n,n1);alpha=zeros(n,n1);

for j=1:n1
    
for i=1:n
    
     Tr(i,j)=T(i)/Tc(j);
        
               alpha(i,j)=(Tr(i,j)^(N(j)*(M(j)-1)))*(exp(L(j)*(1-(Tr(i,j)^(M(j)*N(j))))));
        
        f=@(x)((R*T(i))/(x-b(j)))-((a(j)*alpha(i,j))/(((x+c(j))*(x+b(j)+(2*c(j))))+((b(j)+c(j))*(x-b(j)))))-P;
     
     OF(i,j)=fzero(f,x0(i));
%           OF(i)=fsolve(f,x0(i));
% OF(i)=fminsearch(f,x0(i))

end

end


OBV=zeros(n,n1);
for j=1:n1
    
for i=1:n
    
    OBV(i,j)=OF(i,j)-((R*T(i))/P);
    
end

end

I=zeros(n1,1);
for i=1:n1
    
  I(i)=min(abs(OBV(:,i)));

end

Andis=zeros(n1,1);
for i=1:n1
    
  Andis(i)=find(abs(OBV(:,i))==I(i));

end

T_Final=zeros(n1,1);
for i=1:n1

    T_Final(i)=T(Andis(i));
    
end
