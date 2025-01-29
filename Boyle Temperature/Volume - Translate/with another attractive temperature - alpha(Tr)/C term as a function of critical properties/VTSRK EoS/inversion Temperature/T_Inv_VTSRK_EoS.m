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
     a(i)=(((0.42747*(R^2)*(Tc(i)^2)))/Pc(i));
     b(i)=(0.08664*R*Tc(i))/Pc(i);
    
  c1(i)=(-45.7247*((1/3)-Zc(i)));
  c2(i)=((-2.184*exp(c1(i)))+0.2658);
  c(i)=((((1/3)-Zc(i))*((R*Tc(i))/Pc(i)))*c2(i));
  
end




T_B=[1173.74
1188.12
1287.07
1306.05
1311.23
1366.46
1413.06
1359.08
1398.20
1412.58
1424.09
1432.72
1442.50
1451.70
1470.11
1481.04
1491.40
1502.90
1512.11
1523.04
1532.24
1543.75
1554.10
1563.88
1573.66
1585.17
1594.95
1605.30
1615.08
1624.86
1634.64
1653.63
1663.41
1673.19
1681.82
1691.60
1701.38
1710.01
1101.83
1208.25];

M=input('Enter Name of Alpha function, M=','s');

switch lower(M)
    case 'ozokwelu'
        
        X=zeros(n1,1);
        T_inv=zeros(n1,1);
        k=zeros(n1,1);
        for i=1:n1
  m(i)=0.266+(0.4459*w(i)^0.5);
  n4(i)=(1/m(i))*(0.2469+(0.7495*w(i)));
  
            F=@(x)  -(((((a(i)*(((-2*m(i)*n4(i))/Tc(i))*((x/Tc(i))^(n4(i)-1))*...
                ((exp(2*m(i)*(1-((x/Tc(i))^n4(i))))))))*R*x)-...
                (R*a(i)*(((exp(m(i)*(1-((x/Tc(i))^n4(i)))))^2)))))/...
    ((R*x)^2))-((b(i)-c(i)-((a(i)*(((exp(m(i)*(1-((x/Tc(i))^n4(i)))))^2)))/(R*x)))/x);

T_inv(i)=fzero(F,T_B(i));

X(i)=T_inv(i)/T_B(i);

% T_inv_r(i)=T_inv(i)/Tc(i);
% T_B_r(i)=T_B(i)/Tc(i);
% Pr(i)=P/Pc(i);

        end
        
    case 'nasrifar'
             X=zeros(n1,1);
        T_inv=zeros(n1,1);
        k=zeros(n1,1);
        for i=1:n1
 k(i)=0.480+(1.574*w(i))-(0.176*w(i)^2);
  b1(i)=0.25*(12-(11*k(i))+(k(i)^2));
  b2(i)=0.5*(-6+(9*k(i))-(k(i)^2));
  b3(i)=0.25*(4-(7*k(i))+(k(i)^2));
  
            F=@(x)   -(((((a(i)*(-((Tc(i)*b1(i))/(x^2))-...
        ((2*(Tc(i)^2)*b2(i))/(x^3))-...
        ((3*(Tc(i)^3)*b3(i))/(x^4))))*R*x)-...
                (R*a(i)*((((b1(i)/(x/Tc(i)))+...
         (b2(i)/((x/Tc(i)))^2)+(b3(i)/((x/Tc(i)))^3)))))))/...
    ((R*x)^2))-((b(i)-c(i)-((a(i)*((((b1(i)/(x/Tc(i)))+...
         (b2(i)/((x/Tc(i)))^2)+(b3(i)/((x/Tc(i)))^3)))))/(R*x)))/x);

T_inv(i)=fzero(F,T_B(i));

X(i)=T_inv(i)/T_B(i);

% T_inv_r(i)=T_inv(i)/Tc(i);
% T_B_r(i)=T_B(i)/Tc(i);
% Pr(i)=P/Pc(i);

        end

                

end





