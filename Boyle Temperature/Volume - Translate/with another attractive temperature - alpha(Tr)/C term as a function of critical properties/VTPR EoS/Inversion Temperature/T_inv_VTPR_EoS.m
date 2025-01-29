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
pc=data.Pc;
w=data.w;
Zc=data.Zc;
n1=numel(Tc);


P=1e-2;R=8.3144598;

Pc=pc.*10^5;

b=zeros(n1,1);
a=zeros(n1,1);
for i=1:n1
    
  a(i)=(0.45724*(R^2)*(Tc(i)^2))/Pc(i);
  b(i)=(0.07780*R*Tc(i))/Pc(i);
    c(i)=(-0.252*((R*Tc(i))/Pc(i))*((1.5448*Zc(i))-0.4024));

  
end




T_B=[1005.15
1031.42
1124.91
1157.48
1176.39
1231.54
1279.86
1327.66
1370.21
1395.42
1419.06
1439.54
1461.61
1482.09
1522.01
1543.55
1564.03
1585.04
1603.43
1624.44
1642.30
1663.31
1681.17
1699.55
1717.41
1738.42
1756.80
1775.19
1793.57
1811.96
1829.81
1863.96
1882.87
1901.78
1917.01
1936.44
1955.88
1971.64
973.64
1072.91];

M=input('Enter Name of Alpha function, M=','s');

switch lower(M)
    case 'coquelet'
        
        X=zeros(n1,1);
        T_inv=zeros(n1,1);
        k=zeros(n1,1);
        for i=1:n1
              k(i)=0.41287+(1.34494*w(i))+(0.00421*w(i)^2);

            F=@(x)  -(((((a(i)*(-(k(i)/Tc(i))*(exp(k(i)*(1-(x/Tc(i)))))))*R*x)-...
                (R*a(i)*(exp(k(i)*(1-(x/Tc(i))))))))/...
    ((R*x)^2))-((b(i)-c(i)-((a(i)*exp(k(i)*(1-(x/Tc(i)))))/(R*x)))/x);

T_inv(i)=fzero(F,T_B(i));

X(i)=T_inv(i)/T_B(i);

% T_inv_r(i)=T_inv(i)/Tc(i);
% T_B_r(i)=T_B(i)/Tc(i);
% Pr(i)=P/Pc(i);

        end
        
    case 'original'

                X=zeros(n1,1);
        T_inv=zeros(n1,1);
        k=zeros(n1,1);
        for i=1:n1
              k(i)=0.37464+(1.54226*w(i))-(0.26992*w(i)^2); % original 


            F=@(x)  -(((((a(i)*(-(k(i)/x)*(sqrt(x/Tc(i)))*...
                (sqrt(((1+(k(i)*(1-((x/Tc(i))^0.5))))^2)))))*R*x)-...
                (R*a(i)*((1+(k(i)*(1-((x/Tc(i))^0.5))))^2))))/...
    ((R*x)^2))-((b(i)-c(i)-((a(i)*((1+(k(i)*(1-((x/Tc(i))^0.5))))^2))/(R*x)))/x);

T_inv(i)=fzero(F,T_B(i));

X(i)=T_inv(i)/T_B(i);

% T_inv_r(i)=T_inv(i)/Tc(i);
% T_B_r(i)=T_B(i)/Tc(i);
% Pr(i)=P/Pc(i);

        end

end





