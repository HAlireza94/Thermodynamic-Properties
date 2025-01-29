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

      Oa(i)=(0.66121-(0.76105*Zc(i)));
 Ob(i)=(0.02207+(0.20868*Zc(i)));
   Oc(i)=(0.57765-(1.87080*Zc(i)));

    a(i)=(((Oa(i)*(R^2)*(Tc(i)^2)))/Pc(i));
   b(i)=((Ob(i)*R*Tc(i))/Pc(i));
 c(i)=(Oc(i)*R*Tc(i))/Pc(i);
    

 
end



T_B=[1034.33
964.96
1151.15
1165.77
1181.13
1233.93
1279.85
1326.25
1367.45
1386.42
1404.81
1420.02
1437.02
1453.87
1486.34
1505.92
1525.18
1545.95
1562.23
1583.72
1600.67
1623.22
1641.63
1660.42
1679.12
1703.58
1723.52
1744.07
1764.71
1786
1807.51
1846.4
1869.52
1893.54
1911.02
1935.64
1961.18
1979.76
988.69
1083.18];

                X=zeros(n1,1);
        T_inv=zeros(n1,1);
        
        for i=1:n1    

            k(i)=0.46283+(3.58230*w(i)*Zc(i))+(8.19417*(w(i)^2)*(Zc(i)^2));


            F=@(x)  -(((((a(i)*(-(k(i)/x)*(sqrt(x/Tc(i)))*...
                (sqrt(((1+(k(i)*(1-((x/Tc(i))^0.5))))^2)))))*R*x)-...
                (R*a(i)*((1+(k(i)*(1-((x/Tc(i))^0.5))))^2))))/...
    ((R*x)^2))-((b(i)-((a(i)*((1+(k(i)*(1-((x/Tc(i))^0.5))))^2))/(R*x)))/x);

              T_inv(i)=fzero(F,T_B(i));

              X(i)=T_inv(i)/T_B(i);
              
              Tr_inv(i)=T_inv(i)/Tc(i);
              Tr_B(i)=T_B(i)/Tc(i);

        end







