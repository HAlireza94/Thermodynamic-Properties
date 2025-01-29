clc
clear 
close all



format long

N=input('Enter Group, name=','s');
name=lower(N);
switch name
    case 'alkanols'
        data1=load('Alkanols');
    case 'thiophenes'
        data1=load('Thiophenes');
    case 'pyridines'
        data1=load('Pyridines');
    case 'alkanes'
        data1=load('Alkanes');
        data2=load('Coefficient_Cp');
        case 'polar gases'
        data1=load('Polar gases');
    case 'alkenes'
        data1=load('Alkenes');
        data2=load('Coefficient_Cp_Alkenes');

    case 'cycloalkanes'
        data1=load('Cycloalkanes');
        data2=load('Coefficient_Cp_Cycloalkanes');
    case 'amines'
        data1=load('Amines');
    case 'glycol ethers'
        data1=load('Glycol ethers');
    case 'water'
        data1=load('Water');
    case 'aromatics'
        data1=load('Aromatics');
    case 'gases'
        data1=load('Gases');
    case 'ethers'
        data1=load('Ethers');
    case 'ketones'
        data1=load('Ketones');
        case 'halogens'
        data1=load('Halogenes');
        case 'noble gases'
        data1=load('Noble gases');
        
    case 'fatty acid esters'
        
        data1=load('Fatty acid esters');
     case 'vd'
        data1=load('Validate_data');
        
        data2=load('Coefficient_Cp_Validate');
end


Tc=data1.Tc;
pc=data1.Pc;
w=data1.w;
Zc=data1.Zc;

M=data1.Mw;

Pc=pc.*10^5;

n1=numel(Tc);

P=1e-2;R=8.3144598;

A=data2.A;
B=data2.B;
C=data2.C;
D=data2.D;
E=data2.E;


for i=1:n1
    
  a(i)=(0.45724*(R^2)*(Tc(i)^2))/Pc(i);
  b(i)=(0.07780*R*Tc(i))/Pc(i);
  
  k(i)=0.37464+(1.54226*w(i))-(0.26992*w(i)^2); % original 

end

Pr=P/Pc;



T=[1019.22
1103.80
1222.22
1306.81
1384.88];


 n=numel(T);

 
 x0=zeros(n,1);
for i=1:n
    
    
    x0(i)=(R*T(i))/P;
    
    
end


OF=zeros(n,1);
for i=1:n
    
     f=@(x) ((R*T(i))/(x-b(i)))-((a(i)*((1+(k(i)*(1-((T(i)/Tc(i))^0.5))))^2))/((x^2)+(2*b(i)*x)-(b(i)^2)))-P;
     
     OF(i)=fzero(f,x0(i)); % OF : Volume
%           OF(i)=fsolve(f,x0(i));
% OF(i)=fminsearch(f,x0(i))

end

Cv=zeros(n,1);
Cp=zeros(n,1);
Gamma=zeros(n,1);
X=zeros(n,1);
sound_speed=zeros(n,1);
RSS=zeros(n,1);

for i=1:n
      
    Cp_star(i)=A(i)+(B(i)*T(i))+(C(i)*(T(i)^2))+(D(i)*(T(i)^3))+(E(i)*(T(i)^4));
    
    Cv_star(i)=Cp_star(i)-R;
    
    Tr(i)=T(i)/Tc(i);
    
    Z(i)=(P*OF(i))/(R*T(i));
    
    alpha(i)=((1+(k(i)*(1-((T(i)/Tc(i))^0.5))))^2);
    
    first_order_alpha(i)=-(k(i)/T(i))*(sqrt(Tr(i)))*(sqrt(alpha(i)));
     

      second_order_alpha(i)=(-k(i)*((-k(i)*sqrt(T(i))-...
    (sqrt(Tc(i))*sqrt(alpha(i))))))/(2*Tc(i)*(T(i)^1.5));
   



 Cv(i)=Cv_star(i)+(T(i)*...
     ((a(i)*second_order_alpha(i))/(2*sqrt(2)*b(i)))*(log((OF(i)+...
     ((1+sqrt(2))*b(i)))/(OF(i)+((1-sqrt(2))*b(i))))));
    
    first_order_P_T(i)=(R/(OF(i)-b(i)))-((a(i)*first_order_alpha(i))/((OF(i)^2)+(2*b(i)*OF(i))-(b(i)^2)));
    
    first_order_P_V(i)=(((2*a(i)*alpha(i)*(OF(i)+b(i)))/(((OF(i)^2)+(2*b(i)*OF(i))-(b(i)^2))^2)))-((R*T(i))/((OF(i)-b(i))^2));
    
    Cp(i)=Cv(i)-((T(i)*(first_order_P_T(i)^2))/(first_order_P_V(i)));
    
    X(i)=Cp(i)*Cv(i);
    
    Gamma(i)=Cp(i)/Cv(i);

   SS(i)=sqrt((-OF(i)^2)*(Gamma(i)/M(i))*first_order_P_V(i)); 
    
    sound_speed(i)=SS(i)*sqrt(1000); %$#@
    
    Kappa(i)=((-(1/OF(i)))*(first_order_P_V(i)^(-1)));
    
    
    J_T(i)=((T(i)*((-first_order_P_T(i))/first_order_P_V(i)))-OF(i))/Cp(i);
    
    
end
    
    

figure (1)

plot(T,Cp)
xlabel('Temperature (K)');ylabel('heat capacity at constant pressure(Cp)')

figure (2)

plot(T,Cv)
xlabel('Temperature (K)');ylabel('heat capacity at constant volume (Cv)')


figure (3)
plot(T,Gamma)
xlabel('Temperature (K)');ylabel('Gamma (Cp/Cv)')




