clc
clear 
close all

format long

% Tc=150.86;pc=48.98;w=0.000;Zc=0.291;M=39.948; % Argon

Tc=190.6;pc=46.1;w=0.011;M=16.0; % Methane

% Tc=305.3;Pc=49;w=0.099; % Ethane

% Tc=369.9;Pc=42.5;w=0.153;  % Propane

% Tc=425.12;Pc=37.96;w=0.2;Zc=0.274; % Butane

% Tc=469.7;Pc=33.7;w=0.252; % Pentane

% Tc=507.6;Pc=30.25;w=0.301;Zc=0.266; % Hexane

% Tc=540.2;Pc=27.4;w=0.350; % Heptane

% Tc=568.7;Pc=24.9;w=0.400;Zc=0.256; % Octane

% Tc=594.6;Pc=22.9;w=0.444;Zc=0.235; % Nonane

% Tc=617.7;Pc=21.1;w=0.492; % Decane

%     Tc=788.59;Pc=11.86;w=0.891;Zc=0.210; % Eicosane

%  Tc=914.38;Pc=7.35;w=1.212;Zc=0.176; % teriacontane

%      Tc=1013.50;Pc=5.60;w=1.466;Zc=0.148; % tetracontane

%     Tc=1096.55;Pc=3.95;w=1.678;Zc=0.125; % pentacontane

%     Tc=1168.68;Pc=3.05;w=1.862;Zc=0.111; % hexacontane

%     Tc=1232.84;Pc=2.53;w=2.026;Zc=0.104; % heptacontane

%     Tc=1290.89;Pc=2.13;w=2.174;Zc=0.098; % octacontane

%     Tc=1344.06;Pc=1.81;w=2.310;Zc=0.092; % nonacontane

%     Tc=1393.24;Pc=1.56;w=2.436;Zc=0.086;  % hectane

% Tc=304.2;Pc=73.8;w=0.239;Zc=0.274; % Carbon dioxid

% Tc=562;Pc=48.9;w=0.212;Zc=0.268; % Benzene

% Tc=513;Pc=81;w=0.556;Zc=0.224; % Methanol

% Tc=647.13;Pc=220.55;w=0.34449;Zc=0.229; %Water


% P=1e-7;R=83.14472;

Pc=pc.*10^5;
P=0.01*Pc;
R=8.3144598;


a=(0.45724*(R^2)*(Tc^2))/Pc;b=(0.07780*R*Tc)/Pc;


k=0.41287+(1.34494*w)+(0.00421*w^2); % Coquelet
x1=input('Enter lowe range of temperature, T1=');
x2=input('Enter uper range of temperature, T2=');


 T=linspace(x1,x2,100)';

 n=numel(T);

 
 x0=zeros(n,1);
for i=1:n
    
    
    x0(i)=(R*T(i))/P;
    
    
end


OF=zeros(n,1);
for i=1:n
    
     f=@(x) ((R*T(i))/(x-b))-((a*(exp(k*(1-(T(i)/Tc)))))/((x^2)+(2*b*x)-(b^2)))-P;
     
     OF(i)=fzero(f,x0(i)); % OF : Volume
%           OF(i)=fsolve(f,x0(i));
% OF(i)=fminsearch(f,x0(i))

end

Delta_Helmholtz_Energy=zeros(n,1);
Delta_Internal_Energy=zeros(n,1);
Delta_Entropy=zeros(n,1);Tr=zeros(n,1);
Delta_Enthalpy=zeros(n,1);B=zeros(n,1);
Z=zeros(n,1);F_B=zeros(n,1);alpha=zeros(n,1);
delta_SS=zeros(n,1);delta_Cv=zeros(n,1);
delta_Cp=zeros(n,1);delta_Gamma=zeros(n,1);
delta_sound_speed=zeros(n,1);second_order_alpha=zeros(n,1);
first_order_alpha=zeros(n,1);first_order_P_T=zeros(n,1);
first_order_P_V=zeros(n,1);DEI=zeros(n,1);delta_J_T=zeros(n,1);
for i=1:n
    
    Tr(i)=T(i)/Tc;
    B(i)=(b*P)/(R*T(i));
    Z(i)=(P*OF(i))/(R*T(i));
    
    alpha(i)=exp(k*(1-Tr(i))); % Coquelet
    
    
    
    DEI(i)=(a*alpha(i)-(T(i)*((-k*a*alpha(i))/Tc)));
    
    
    Delta_Enthalpy(i)=((((DEI(i))/(2*sqrt(2)*b)))*...
        (log((OF(i)+((1-sqrt(2))*b))/(OF(i)+((1+sqrt(2))*b)))))+(R*T(i)*(Z(i)-1)); % Coquelet
    
    Delta_Entropy(i)=(R*(log(Z(i)-B(i))))+...
        (((((k*a)/Tc)*alpha(i))/(2*sqrt(2)*b))*...
        (log((OF(i)+((1-sqrt(2))*b))/(OF(i)+((1+sqrt(2))*b))))); % Coquelet
    
    
    Delta_Helmholtz_Energy(i)=(R*T(i)*(log(1/(Z(i)-B(i)))))+...
        (((a*alpha(i))/(2*sqrt(2)*b))*...
        (log((OF(i)+((1-sqrt(2))*b))/(OF(i)+((1+sqrt(2))*b)))));
    
    
    Delta_Internal_Energy(i)=((((DEI(i))/(2*sqrt(2)*b)))*...
        (log((OF(i)+((1-sqrt(2))*b))/(OF(i)+((1+sqrt(2))*b)))));

    
    F_B(i)=Delta_Enthalpy(i)/Delta_Entropy(i);
   
    
 
    
end




I=max(F_B);
disp('Objective Value, alpha=');disp(I);

Andis=find(F_B==I);
% [row,colm]=find(OBV==I)

T_Final=T(Andis);
disp('Temperature, T');disp(T_Final)

RH=Delta_Enthalpy(Andis);
RS=Delta_Entropy(Andis);


plot(T_Final,I,'ok',T,F_B)

legend({'Boyle Temperature','Ratio of Residual Enthalpy to Residual Entropy'})
