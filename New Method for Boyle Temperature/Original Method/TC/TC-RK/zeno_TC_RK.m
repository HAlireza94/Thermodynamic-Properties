clc
clear
close all
format long

% Tc=150.86;Pc=48.98;w=0;Zra=0.2941;L=0.202;M=0.9085;N=1.8138;Vc=74.59;TB=351.2631316; %Ar
% Tc=304.21;Pc=73.83;w=0.2236;Zra=0.2728;L=0.2806;M=0.8684;N=2.2782;Vc=94.0;TB=619.805903; %CO2
% Tc=126.2;Pc=34;w=0.0377;Zra=0.2908;L=0.1901;M=0.89;N=2.0107;Vc=89.21;TB=285.3591796; %N2
Tc=190.564;Pc=45.99;w=0.0115;Zra=0.2894;L=0.217;M=0.9082;N=1.8172;Vc=98.6;TB=442.2161081; %methane


R=83.144598;Z=1;


a=(1/(9*((2^(1/3))-1)))*((((R^2)*(Tc^2)))/Pc);
   

            c=((R*Tc)/Pc)*(0.2150-(0.7314*Zra));
        
        b=((((2^(1/3))-1)/3)*((R*Tc)/Pc))-c;


  RhoR=linspace(1e-10,3,100)';
%                  x0=linspace(2,0.5,100)';
                 x0=linspace((TB/Tc),0.5,100);
                 Tr=zeros(100,1);
                 for i=1:100
                             


A(i)=((R*Tc*RhoR(i))/(Vc-(b*RhoR(i))));
B(i)=((Z*R*Tc*RhoR(i))/Vc);
C(i)=((a)/((((Vc/RhoR(i))+c)*((Vc/RhoR(i))+b+(2*c)))));
E(i)=C(i)*exp(L);
D=N*(M-1);

% f=@(x)((A(i)-B(i))*x)-(C(i)*((x^(N*(M-1)))*(exp(L*(1-(x^(M*N)))))));

f=@(x)((A(i)-B(i))*x)-(E(i)*(x^D)*exp(-L*(x^(M*N))));

Tr(i)=fzero(f,x0(i));

                  

                 end
                 
                 plot(RhoR,Tr,'+');

    KM=[RhoR,Tr];

