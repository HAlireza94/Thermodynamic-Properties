clc
clear
close all
format long

% Zra=0.2941;L=0.1227;M=0.9045;N=1.8541;Tc=150.86;Pc=48.98;w=0;Vc=74.59;TB=383.6787679; %Argon
% Tc=304.21;Pc=73.83;w=0.2236;Zra=0.2728;L=0.1784;M=0.859;N=2.4107;Vc=94.0;TB=659.1115112; % CO2
% Tc=126.2;Pc=34;w=0.0377;Zra=0.2908;L=0.1242;M=0.8898;N=2.013;Vc=89.21;TB=3.09E+02; %N2
Tc=190.564;Pc=45.99;w=0.0115;Zra=0.2894;L=0.1473;M=0.9075;N=1.8243;Vc=98.6;TB=477.7289729; %methane


R=83.144598;Z=1;


a=(0.45724*(R^2)*(Tc^2))/Pc;
c=((R*Tc)/Pc)*(0.1975-(0.7325*Zra));
        b=((0.07780*R*Tc)/Pc)-c;

  RhoR=linspace(1e-10,3,100)';
%                  x0=linspace(2,0.5,100)';
                 x0=linspace((TB/Tc),0.5,100);
                 Tr=zeros(100,1);
                 for i=1:100
                             

%                     f=@(x) (((R*x*Tc*RhoR(i))/(Vc-(b*RhoR(i))))...
%                         -((a*((x^(N*(M-1)))*(exp(L*(1-(x^(M*N)))))))/...
%                         ((((Vc/RhoR(i))+c)*((Vc/RhoR(i))+b+(2*c)))+((b+c)*((Vc/RhoR(i))-b))))...
%                         -((Z*R*x*Tc*RhoR(i))/Vc));

A(i)=((R*Tc*RhoR(i))/(Vc-(b*RhoR(i))));
B(i)=((Z*R*Tc*RhoR(i))/Vc);
C(i)=((a)/((((Vc/RhoR(i))+c)*((Vc/RhoR(i))+b+(2*c)))+((b+c)*((Vc/RhoR(i))-b))));
E(i)=C(i)*exp(L);
D=N*(M-1);

% f=@(x)((A(i)-B(i))*x)-(C(i)*((x^(N*(M-1)))*(exp(L*(1-(x^(M*N)))))));

f=@(x)((A(i)-B(i))*x)-(E(i)*(x^D)*exp(-L*(x^(M*N))));

Tr(i)=fzero(f,x0(i));

                  

                 end
                 
                 plot(RhoR,Tr,'+');

KM=[RhoR,Tr];
