clc
clear
close all
format long 
M=input('Enter Name Of Gas :','s');
p=180;

[Tc, Pc, w]=Critical_Properties(M);
Experimental_Data=ED(M,p);

% A=34.942;B=-0.039957;C=0.00019184;D=-1.5303e-07;E=3.9321e-11;Mw=16;x1=100;x2=550; % Methane

% A=28.146;B=0.043447;C=0.00018946;D=-1.9082e-07;E=5.3349e-11;Mw=30.07;x1=150;x2=600; % Ethane

% A=31.78;B=0.74489;C=-0.00010945;D=-2.2668e-07;E=9.3458e-11;Mw=142.285;x1=500;x2=1000; % Decane

% A=32.083;B=-0.014831;C=0.00024774;D=-2.3766e-07;E=6.8274e-11;Mw=28.054;x1=150;x2=550; % Ethene

% A=31.298;B=0.072449;C=0.00019481;D=-2.1582e-07;E=6.2974e-11;Mw=42.081;x1=150;x2=600; % Propene

% A=29.414;B=-0.0045993;C=0.000013004;D=-5.4759e-09;E=2.9239e-13;Mw=28.013;x1=100;x2=550; % Nitrogen

% A=29.526;B=-8.8999e-03;C=3.8083e-05;D=-3.2629e-08;E=8.8607e-12;Mw=31.999;x1=100;x2=550; % Oxygen
 
% A=27.437;B=0.042315;C=-0.000019555;D=3.9968e-09;E=-2.9872e-13;Mw=44.01;x1=200;x2=550; % CO2

% A=29.556;B=-0.0065807;C=0.00002013;D=-1.2227e-08;E=2.2617e-12;Mw=28.01;x1=100;x2=550; % CO

%  A=20.786;B=0;C=0;D=0;E=0;Mw=39.948;x1=100;x2=550; % Argon

%  A=20.786;B=0;C=0;D=0;E=0;Mw=131.29;x1=150;x2=600; % Xenon

% A=15.278;B=0.19916;C=-0.00016369;D=5.1686e-08;E=-3.182e-12;Mw=88.005;x1=100;x2=550; % Tetrafluoromethane ( R - 14 )

% A=45.579;B=0.49467;C=-0.00040808;D=1.3789e-07;E=-1.1769e-11;Mw=200.032;x1=250;x2=700; % Octafluorocyclobutane (RC-318)

A=-31.368;B=4.75e-01;C=-3.11e-04;D=8.52e-08;E=-5.05e-12;Mw=78.114;x1=400;x2=850; % Benzene

% A=40.046;B=-0.038287;C=0.00024529;D=-2.1679E-07;E=5.9909E-11;Mw=32.042;x1=400;x2=850; % Methanol

% A=33.573;B=-0.012581;C=0.000088906;D=-7.1783e-08;E=1.8569e-11;Mw=17.031;x1=360;x2=810; % Ammonia

% A=33.933;B=-0.0084186;C=0.000029906;D=-1.7825e-08;E=3.6934e-12;Mw=18.015;x1=500;x2=950; % Water


disp('==============================================');

alpha_function='jub';

disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');

jtx=Experimental_Data(:,1);Tx=Experimental_Data(:,2);

% T=linspace(x1,x2,300)'; % K
T=Tx;

P=180*101325;
a0=input('Enter intial number, a0=');
j=numel(T);


s1=numel(P);
J_T=zeros(j,s1);
Cp=zeros(j,s1);
Cv=zeros(j,s1);
Speed_Sound=zeros(j,s1);

for s=1:s1
    
for i=1:j
    
   [J_T(i,s) ,Speed_Sound(i,s), Cp(i,s), Cv(i,s)]=Caclulation_Joule_Thomson_Coefficient_PR_EoS_For_Pure...
       (alpha_function,Tc,Pc,w,T(i),P(s),A,B,C,D,E,a0,Mw);

end

end

jt=J_T.*101325;
plot(T,jt,Tx,jtx,'*')


Err=(jt-jtx).^2;

Answer=(mean(Err))^0.5