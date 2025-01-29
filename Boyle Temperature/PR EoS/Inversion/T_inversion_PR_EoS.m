clc
clear
close all

format long

% N=input('Enter Group, name=','s');
% name=lower(N);
% switch name
%     case 'alkanols'
%         data=load('Alkanols');
%     case 'thiophenes'
%         data=load('Thiophenes');
%     case 'pyridines'
%         data=load('Pyridines');
%     case 'alkanes'
%         data=load('Alkanes');
%     case 'alkenes'
%         data=load('Alkenes');
%     case 'cycloalkanes'
%         data=load('Cycloalkanes');
%     case 'amines'
%         data=load('Amines');
%     case 'glycol ethers'
%         data=load('Glycol ethers');
%     case 'water'
%         data=load('Water');
%     case 'aromatics'
%         data=load('Aromatics');
%     case 'gases'
%         data=load('Gases');
%     case 'ethers'
%         data=load('Ethers');
%     case 'ketones'
%         data=load('Ketones');
%          case 'halogenes'
%         data=load('Halogenes');
%         case 'noble gases'
%         data=load('Noble gases');
%         case 'polar gases'
%         data=load('Polar gases');
% end
% 
% Tc=data.Tc;
% pc=data.Pc;
% w=data.w;

Tc=[430.05
456.15
496.95
531.9
557.66
589.15
617.7
643.92
668.22
690.9
712.21
732.33
751.41
769.57
803.53
819.48
834.82
849.62
863.91
877.74
891.14
904.15
916.78
929.06
941.02
952.68
964.05
975.16
986.01
996.61
1006.99
1027.11
1036.88
1046.45
1055.85
1065.07
1074.13
1083.04];

pc=[74.6
56.2
47.4
42
37.28
32.88
29.49
26.61
24.2
22.18
20.41
18.95
17.59
16.39
14.36
14.21
13.54
11.72
10.98
11.03
10.49
9.99
9.52
9.07
8.65
8.25
7.88
7.53
7.63
7.32
7.02
6.48
6.23
5.99
5.77
5.55
5.35
5.16];


w=[0.281
0.285
0.28
0.329
0.364
0.399
0.431
0.46
0.487
0.538
0.592
0.644
0.692
0.739
0.826
0.866
0.906
0.943
0.98
1.015
1.05
1.083
1.115
1.146
1.177
1.207
1.236
1.264
1.292
1.319
1.346
1.397
1.422
1.446
1.47
1.494
1.517
1.54];

n1=numel(Tc);


P=1e-2;R=8.3144598;

Pc=pc.*10^5;

b=zeros(n1,1);
a=zeros(n1,1);
for i=1:n1
    
  a(i)=(0.45724*(R^2)*(Tc(i)^2))/Pc(i);
  b(i)=(0.07780*R*Tc(i))/Pc(i);
  
end




T_B=[990.29
1048.65
1145.55
1197.30
1234.73
1285.39
1330.53
1371.27
1408.71
1430.73
1449.45
1467.07
1484.68
1501.20
1534.23
1550.75
1565.07
1580.48
1594.79
1610.21
1623.42
1637.74
1650.95
1664.16
1676.28
1689.49
1701.60
1713.71
1725.83
1736.84
1747.85
1770.97
1780.88
1791.89
1802.90
1812.81
1822.72
1832.63];

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
    ((R*x)^2))-((b(i)-((a(i)*exp(k(i)*(1-(x/Tc(i)))))/(R*x)))/x);

T_inv(i)=fzero(F,T_B(i));

X(i)=T_inv(i)/T_B(i);

Tr_inv(i)=T_inv(i)/Tc(i);
Tr_B(i)=T_B(i)/Tc(i);
Pr(i)=P/Pc(i);


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
    ((R*x)^2))-((b(i)-((a(i)*((1+(k(i)*(1-((x/Tc(i))^0.5))))^2))/(R*x)))/x);

T_inv(i)=fzero(F,T_B(i));

X(i)=T_inv(i)/T_B(i);

Tr_inv(i)=T_inv(i)/Tc(i);
Tr_B(i)=T_B(i)/Tc(i);
Pr(i)=P/Pc(i);

        end
        
    case 'haghtalab'
        
                X=zeros(n1,1);
        T_inv=zeros(n1,1);
        n2=zeros(n1,1);
        for i=1:n1
             n2(i)=1.7309+(1.6571*w(i))+(0.1554*w(i)^2); % Haghtalab et al. (2010)

            F=@(x)  -(((((a(i)*((-(exp(1-(n2(i)^(log(x/Tc(i))))))*(log(n2(i)))*((x/Tc(i))^(log(n2(i)))))/x))*R*x)-...
                (R*a(i)*(exp(1-(n2(i)^(log(x/Tc(i)))))))))/...
                 ((R*x)^2))-((b(i)-((a(i)*(exp(1-(n2(i)^(log(x/Tc(i)))))))/(R*x)))/x);
T_inv(i)=fzero(F,T_B(i));

X(i)=T_inv(i)/T_B(i);

Tr_inv(i)=T_inv(i)/Tc(i);
Tr_B(i)=T_B(i)/Tc(i);
Pr(i)=P/Pc(i);


        end


end





