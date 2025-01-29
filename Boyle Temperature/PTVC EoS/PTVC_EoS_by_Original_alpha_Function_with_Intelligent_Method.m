clc
clear 
close all

format long

N=input('Enter Group, name=','s');
name=lower(N);
switch name
    case 'alkanols'
        data=load('Alkanols');
    case 'thiophenes'
        data=load('Thiophenes');
    case 'pyridines'
        data=load('Pyridines');
    case 'alkanes'
        data=load('Alkanes');
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
    case 'ehters'
        data=load('Ethers');
    case 'ketones'
        data=load('Ketones');
         case 'halogenes'
        data=load('Halogenes');
    case 'noble gases'
        data=load('Noble gases');
    case 'polar gases'
        data=load('Polar gases');
end


Tc=data.Tc;
Pc=data.Pc;
w=data.w;
Zc=data.Zc;

n1=numel(Tc);
P=1e-07;R=83.14472;
F=zeros(n1,1);
c=zeros(n1,1);
b=zeros(n1,1);
a=zeros(n1,1);
Oc=zeros(n1,1);
Ob=zeros(n1,1);
Oa=zeros(n1,1);

for i=1:n1
  
    Oa(i)=0.69368018-(1.0634424*Zc(i))+( 0.68289995*Zc(i)^2)+(-0.21044403*Zc(i)^3)+(0.003752658*Zc(i)^4);
 Ob(i)=0.025987178+(0.180754784*Zc(i))+(0.061258949*Zc(i)^2);
   Oc(i)=0.577500514-(1.898414283*Zc(i));

    a(i)=(((Oa(i)*(R^2)*(Tc(i)^2)))/Pc(i));
   b(i)=((Ob(i)*R*Tc(i))/Pc(i));
 c(i)=(Oc(i)*R*Tc(i))/Pc(i);
    F(i)=-6.608+(70.43*Zc(i))-(159.0*Zc(i)^2);  

  
end

x1=input('Enter lowe range of temperature, T1=');
x2=input('Enter uper range of temperature, T2=');


 T=linspace(x1,x2,100000)';n=numel(T);

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


