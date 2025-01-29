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
    case 'ethers'
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
P=1e-07;
n1=numel(Tc);
R=83.14472;
n3=zeros(n1,1);
b=zeros(n1,1);
a=zeros(n1,1);
for i=1:n1
    
  a(i)=(0.45724*(R^2)*(Tc(i)^2))/Pc(i);
  b(i)=(0.07780*R*Tc(i))/Pc(i);
  n3(i)=1.7309+(1.6571*w(i))+(0.1554*w(i)^2);


end




x1=input('Enter lowe range of temperature, T1=');
x2=input('Enter uper range of temperature, T2=');


 T=linspace(x1,x2,1000)';n=numel(T);

 
 x0=zeros(n,1);
for i=1:n
    
    x0(i)=(R*T(i))/P;
    
end


OF=zeros(n,n1);
for j=1:n1
    
for i=1:n
    Tr(i,j)=T(i)/Tc(j);
    B(i,j)=(b(j)*P)/(R*T(i)); 
    alpha(i,j)=(exp(1-(n3(j)^(log(Tr(i,j)))))); % Haghtalab et al. (2011)
    f=@(x) ((R*T(i))/(x-b(j)))-((a(j)*(exp(1-(n3(j)^(log(T(i)/Tc(j)))))))/((x^2)+(2*b(j)*x)-(b(j)^2)))-P;
     
     OF(i,j)=fzero(f,x0(i));
%           OF(i)=fsolve(f,x0(i));
% OF(i)=fminsearch(f,x0(i))

end

end


for i=1:n
    
    Z(i)=(P*OF(i))/(R*T(i));
   
end


CS=zeros(n,n1);
Delta_Enthalpy=zeros(n,n1);
Delta_Entropy=zeros(n,n1);
for j=1:n1
    
for i=1:n
  
    Delta_Enthalpy(i,j)=(R*T(i)*(Z(i)-1))+...
        ((((a(j)*alpha(i,j))+(a(j)*(n3(j)^log(Tr(i,j))*alpha(i,j)*log(n3(j)))))/(2*sqrt(2)*b(j)))...
        *(log((OF(i)+((1-sqrt(2))*b(j)))/(OF(i)+((1+sqrt(2))*b(j))))));
    
    
    
    
    Delta_Entropy(i,j)=(R*(log(Z(i)-B(i,j))))+...
        ((((n3(j)^log(Tr(i,j))*a(j)*alpha(i,j)*log(n3(j)))/T(i))/(2*sqrt(2)*b(j)))*...
        (log((OF(i)+((1-sqrt(2))*b(j)))/(OF(i)+((1+sqrt(2))*b(j))))));

    cs(i,j)=((1/log(n3(j)))-1)*((((n3(j)^log(Tr(i,j))*a(j)*alpha(i,j)*log(n3(j)))/T(i))/(2*sqrt(2)*b(j)))*...
        (log((OF(i)+((1-sqrt(2))*b(j)))/(OF(i)+((1+sqrt(2))*b(j))))));
    
    
    
    CS(i,j)=(Delta_Entropy(i,j)+(cs(i,j)));
   
    fx(i,j)=Delta_Enthalpy(i,j)/CS(i,j);
    

    
end
    
end
    
    I=zeros(n1,1);
for i=1:n1
    
  I(i)=max(fx(:,i));

end

% disp('Objective Value, alpha=');disp(I);



Andis=zeros(n1,1);
for i=1:n1
    
  Andis(i)=find((fx(:,i))==I(i));

end




T_Final=zeros(n1,1);
for i=1:n1

    T_Final(i)=T(Andis(i));
    
end

% disp('Temperature, T');disp(T_Final)



% plot(T_Final,I,'+',T,fx)

    