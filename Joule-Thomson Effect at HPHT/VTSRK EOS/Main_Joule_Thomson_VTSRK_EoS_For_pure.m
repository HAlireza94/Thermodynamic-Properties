clc
clear
close all

M=input('Enter Name Of Gas :','s');

[Tc, Pc, w, Zc]=Critical_Properties(M);


% A=34.942;B=-0.039957;C=0.00019184;D=-1.5303e-07;E=3.9321e-11;Mw=16; % Methane

% A=28.146;B=0.043447;C=0.00018946;D=-1.9082e-07;E=5.3349e-11; % Ethane

% A=31.78;B=0.74489;C=-0.00010945;D=-2.2668e-07;E=9.3458e-11; % Decane

% A=32.083;B=-0.014831;C=0.00024774;D=-2.3766e-07;E=6.8274e-11; % Ethene

% A=31.298;B=0.072449;C=0.00019481;D=-2.1582e-07;E=6.2974e-11; % Propene

% A=29.414;B=-0.0045993;C=0.000013004;D=-5.4759e-09;E=2.9239e-13; % Nitrogen

% A=29.526;B=-8.8999e-03;C=3.8083e-05;D=-3.2629e-08;E=8.8607e-12; % Oxygen
 
A=27.437;B=0.042315;C=-0.000019555;D=3.9968e-09;E=-2.9872e-13; % CO2

% A=29.556;B=-0.0065807;C=0.00002013;D=-1.2227e-08;E=2.2617e-12; % CO

%  A=20.786;B=0;C=0;D=0;E=0; % Argon

%  A=20.786;B=0;C=0;D=0;E=0; % Xenon

% A=15.278;B=0.19916;C=-0.00016369;D=5.1686e-08;E=-3.182e-12; % Tetrafluoromethane ( R - 14 )

% A=45.579;B=0.49467;C=-0.00040808;D=1.3789e-07;E=-1.1769e-11; % Octafluorocyclobutane (RC-318)

% A=-31.368;B=4.75e-01;C=-3.11e-04;D=8.52e-08;E=-5.05e-12; % Benzene

% A=40.046;B=-0.038287;C=0.00024529;D=-2.1679E-07;E=5.9909E-11; % Methanol

% A=33.573;B=-0.012581;C=0.000088906;D=-7.1783e-08;E=1.8569e-11; % Ammonia

% A=33.933;B=-0.0084186;C=0.000029906;D=-1.7825e-08;E=3.6934e-12; % Water

disp('==============================================')

% alpha_function=input('name of Alpha Function, alpha :','s');

disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx')

a0=input('Enter intial number, a0=');

T=linspace(200,550,300)';

% T=[];

j=numel(T);

% P=[300*101325 360*101325];
P=300*101325;
s1=numel(P);

J_T=zeros(j,5,s1);
Cp=zeros(j,5,s1);
Cv=zeros(j,5,s1);
Speed_Sound=zeros(j,5,s1);


for q=1:s1
    
for i=1:5
    if i==1
    
        alpha_function='original';
    
    elseif i==2
        
        alpha_function='ozokwelu';
        
    elseif i==3
        
        alpha_function='trebble';
        
    elseif i==4
        
        alpha_function='nasrifar';
        
    elseif i==5
        
        alpha_function='jub';

    end
    
    for g=1:j
    
   [J_T(g,i,q), Speed_Sound(g,i,q), Cp(g,i,q), Cv(g,i,q)]=...
       Caclulation_Joule_Thomson_Coefficient_VTSRK_EoS_For_Pure...
       (alpha_function,Tc,Pc,w,T(g),P(q),A,B,C,D,E,Zc,a0,Mw);
   
   
    end

    
end

end

jt=J_T.*101325;

if numel(T)==300 || numel(T)==1000
    figure (1)
    plot(T,jt(:,1,1),T,jt(:,2,1),T,jt(:,3,1),T,jt(:,4,1),T,jt(:,5,1))
    legend({'Original','ozokwelu','trebble','nasrifar','Jaubert'});title('P=100 atm')
    figure (2)
    plot(T,Cp(:,1,1),T,Cp(:,2,1),T,Cp(:,3,1),T,Cp(:,4,1),T,Cp(:,5,1))
    legend({'Original','ozokwelu','trebble','nasrifar','Jaubert'});title('P=100 atm')
    figure (3)
    plot(T,Cv(:,1,1),T,Cv(:,2,1),T,Cv(:,3,1),T,Cv(:,4,1),T,Cv(:,5,1))
    legend({'Original','ozokwelu','trebble','nasrifar','Jaubert'});title('P=100 atm')
    figure (4)
    plot(T,Speed_Sound(:,1,1),T,Speed_Sound(:,2,1),T,Speed_Sound(:,3,1),T,Speed_Sound(:,4,1),T,Speed_Sound(:,5,1))
    legend({'Original','ozokwelu','trebble','nasrifar','Jaubert'});title('P=100 atm')
    
    
%     figure (2)
%     plot(T,jt(:,1,2),T,jt(:,2,2),T,jt(:,3,2),T,jt(:,4,2))
%     legend({'Original','ozokwelu','trebble','nasrifar'});title('P=360 atm')
    
else 
    
    figure (1)
    plot(T,jt(:,1,1),'k+',T,jt(:,2,1),'g*',T,jt(:,3,1),'bo',T,jt(:,4,1),'r>')
    legend({'Original','ozokwelu','trebble','nasrifar'});title('P=300 atm')
    
    figure (2)
    plot(T,jt(:,1,2),'k+',T,jt(:,2,2),'g*',T,jt(:,3,2),'bo',T,jt(:,4,2),'r>')
    legend({'Original','ozokwelu','trebble','nasrifar'});title('P=360 atm')
        
end

M360=jt(:,:,2);
M300=jt(:,:,1);