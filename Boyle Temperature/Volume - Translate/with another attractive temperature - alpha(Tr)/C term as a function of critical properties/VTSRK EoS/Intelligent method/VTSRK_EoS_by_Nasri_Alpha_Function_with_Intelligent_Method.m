clc
clear 
close all

format long

N=input('Enter Group, name=','s');
name=lower(N);
switch name
    case 'alkanols'
        data=load('Alkanols');
         case 'polar gases'
        data=load('Polar gases');
    case 'thiophenes'
        data=load('Thiophenes');
    case 'pyridines'
        data=load('Pyridines');
    case 'alkanes'
        data=load('Alkanesa');
    case 'alkenes'
        data=load('Alkenes');
    case 'cycloalkanes'
        data=load('Cycloalkanes');
    case 'amines'
        data=load('Amines');
    case 'glycol ethers'
        data=load('Glycol ethers');
    case 'water'
        data=load('Water');
    case 'aromatics'
        data=load('Aromatics');
    case 'gases'
        data=load('Gases');
    case 'ethers'
        data=load('Ethers');
    case 'ketones'
        data=load('Ketones');
        case 'halogens'
        data=load('Halogenes');
        case 'noble gases'
        data=load('Noble gases');
end

Tc=data.Tc;
Pc=data.Pc;
w=data.w;
Zc=data.Zc;


n1=numel(Tc);
P=1e-07;R=83.14472;

b3=zeros(n1,1);
b2=zeros(n1,1);
b1=zeros(n1,1);
k=zeros(n1,1);

c=zeros(n1,1);
c2=zeros(n1,1);
c1=zeros(n1,1);

b=zeros(n1,1);
a=zeros(n1,1);
for i=1:n1
    
  a(i)=(((0.42747*(R^2)*(Tc(i)^2)))/Pc(i));
  b(i)=(0.08664*R*Tc(i))/Pc(i);
  
  c1(i)=(-45.7247*((1/3)-Zc(i)));
  c2(i)=((-2.184*exp(c1(i)))+0.2658);
  c(i)=((((1/3)-Zc(i))*((R*Tc(i))/Pc(i)))*c2(i));
  
  k(i)=0.480+(1.574*w(i))-(0.176*w(i)^2);
  b1(i)=0.25*(12-(11*k(i))+(k(i)^2));
  b2(i)=0.5*(-6+(9*k(i))-(k(i)^2));
  b3(i)=0.25*(4-(7*k(i))+(k(i)^2));

end

x1=input('Enter lowe range of temperature, T1=');
x2=input('Enter uper range of temperature, T2=');


 T=linspace(x1,x2,2000)';n=numel(T);

 x0=zeros(n,1);
for i=1:n
    
    
    x0(i)=(R*T(i))/P;
    
    
end


OF=zeros(n,n1);
for j=1:n1
    
for i=1:n
    
     f=@(x) ((R*T(i))/(x+c(j)-b(j)))-((a(j)*((b1(j)/(T(i)/Tc(j)))+(b2(j)/((T(i)/Tc(j)))^2)+(b3(j)/((T(i)/Tc(j)))^3)))/((x+c(j))*(x+c(j)+b(j))))-P;
     
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
