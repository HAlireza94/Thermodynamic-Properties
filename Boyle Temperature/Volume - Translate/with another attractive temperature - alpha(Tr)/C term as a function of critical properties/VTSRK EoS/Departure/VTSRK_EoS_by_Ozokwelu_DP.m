clc
clear 
close all
format long

% Tc=150.86;Pc=48.98;w=0.000;Zc=0.291; % Argon

% Tc=190.56;Pc=45.99200;w=0.011;Zc=0.2863; %Methane

% Tc=305.3;Pc=49;w=0.099;Zc=0.279; % Ethane

% Tc=369.9;Pc=42.5;w=0.153;Zc=0.276; % Propane

% Tc=425.12;Pc=37.96;w=0.2;Zc=0.274; % Butane

% Tc=469.7;Pc=33.7;w=0.252;Zc=0.270; % Pentane

% Tc=507.6;Pc=30.25;w=0.301;Zc=0.266; % Hexane

% Tc=540.2;Pc=27.4;w=0.350;Zc=0.261; % Heptane1

% Tc=568.7;Pc=24.9;w=0.400;Zc=0.256; % Octane

% Tc=594.6;Pc=22.9;w=0.444;Zc=0.235; % Nonane1

%   Tc=617.7;Pc=21.1;w=0.492;Zc=0.231; % Decane

% Tc=304.2;Pc=73.8;w=0.239;Zc=0.274; % Carbon dioxid

% Tc=562;Pc=48.9;w=0.212;Zc=0.268; % Benzene

% Tc=513;Pc=81;w=0.556;Zc=0.224; % Methanol

% Tc=647.13;Pc=220.55;w=0.34449;Zc=0.229; %Water

n10=input('Enter name :','s');


switch n10
     case 'Ethene'
    Tc=282.5;Pc=50.42;w=0.087;Zc=0.281;
    case 'Propene'
    Tc=365;Pc=46.1;w=0.144;Zc=0.281;  
    case '1-Butene'
    Tc=419.59;Pc=40.2;w=0.187;Zc=0.276;
    case '1-Pentene'
    Tc=464.8;Pc=35.6;w=0.236;Zc=0.275;
    case '1-Hexene' %C2H5OH
    Tc=504;Pc=32.1;w=0.285;Zc=0.272;
    case '1-Heptene'
    Tc=537.30;Pc=29.20;w=0.344;Zc=0.267;
    case '1-Octene'
    Tc=567.00;Pc=26.80;w=0.392;Zc=0.266;
    case '1-Nonene'
    Tc=593.48;Pc=23.30;w=0.374;Zc=0.248;
    case '1-Decene'
    Tc=618.16;Pc=22.20;w=0.417;Zc=0.252;
    case '1-Eicosene'
    Tc=801.49;Pc=11.07;w=0.885;Zc=0.211;
    case '1-Triacontene'
    Tc=927.80;Pc=6.88;w=1.208;Zc=0.177;
    case '1-Tetracontene'
    Tc=1027.21;Pc=5.23;w=1.462;Zc=0.149;
    case '1-Pentacontene'
    Tc=1110.43;Pc=3.69;w=1.675;Zc=0.125;
    case '1-Hexacontene'
    Tc=1182.70;Pc=2.84;w=1.860;Zc=0.111;
    case '1-Heptacontene'
    Tc=1246.96;Pc=2.35;w=2.024;Zc=0.104;
    case '1-Octacontene'
    Tc=1305.08;Pc=1.98;w=2.172;Zc=0.098;
end



P=1e-07;
R=83.14472;
a=(((0.42747*(R^2)*(Tc^2)))/Pc);
b=((0.08664*R*Tc)/Pc);


c1=(-45.7247*((1/3)-Zc));
c2=((-2.184*exp(c1))+0.2658);
c=((((1/3)-Zc)*((R*Tc)/Pc))*c2);

m=0.266+(0.4459*w^0.5);
S=(1/m)*(0.2469+(0.7495*w));


% x1=input('Enter lowe range of temperature, T1=');
x2=input('Enter uper range of temperature, T2=');


 T=linspace(Tc,x2,100)';n=numel(T);

 x0=zeros(n,1);
for i=1:n
    
    
    x0(i)=(R*T(i))/P;
    
    
end


OF=zeros(n,1);
for i=1:n
    
     f=@(x) ((R*T(i))/(x+c-b))-((a*((exp(m*(1-((T(i)/Tc)^S))))^2))/((x+c)*(x+c+b)))-P;
     
     OF(i)=fzero(f,x0(i));
%           OF(i)=fsolve(f,x0(i));
% OF(i)=fminsearch(f,x0(i))

end

OBV=zeros(n,1);
for i=1:n
    
    OBV(i)=OF(i)-((R*T(i))/P);
    
end


Delta_Helmholtz_Energy=zeros(n,1);
Delta_Internal_Energy=zeros(n,1);
Delta_Entropy=zeros(n,1);
Delta_Enthalpy=zeros(n,1);
Z=zeros(n,1);
F_B=zeros(n,1);

for i=1:n
    
    Tr(i)=T(i)/Tc;
    B(i)=(b*P)/(R*T(i));
    C(i)=(c*P)/(R*T(i));
    Z(i)=(P*OF(i))/(R*T(i));
    
    alpha(i)=((exp(m*(1-((T(i)/Tc)^S))))^2);
    
    A(i)=(a*alpha(i))/((R*T(i))^2);
    
    
    
    DEI(i)=(a*alpha(i)-(T(i)*(((-2*m*S*a)/Tc)*(Tr(i)^(S-1))*exp(2*m*(1-(Tr(i)^S))))));
    
    
    Delta_Enthalpy(i)=(R*T(i)*(Z(i)-1))-...
        (((DEI(i))/(b))...
        *(log((b/(OF(i)+c))+1)));
    
    Delta_Entropy(i)=(R*(log(Z(i)+C(i)-B(i))))+...
        ((((((-2*m*S*a)/Tc)*(Tr(i)^(S-1))*exp(2*m*(1-(Tr(i)^S)))))/(b))...
        *(log((b/(OF(i)+c))+1)));
    
    Delta_Internal_Energy(i)=Delta_Enthalpy(i)-(R*T(i)*(Z(i)-1));
    
   Delta_Helmholtz_Energy(i)=((R*T(i))*(log(1/(Z(i)+C(i)-B(i)))))...
        +(((a*alpha(i))/b)*(log((b/(OF(i)+c))+1)));
    
    
    F_B(i)=Delta_Enthalpy(i)/Delta_Entropy(i);
    
    Delta_Gibbps(i)=Delta_Enthalpy(i)-(T(i)*Delta_Entropy(i));
    
end


% FB=zeros(n,1);
% for i=1:n
%     if F_B(i)>0
%         FB(i)=F_B(i);
%     else 
%         FB(i)=0;
%     end
% end

I=max(F_B);
disp('Objective Value, alpha=');disp(I);

Andis=find(F_B==I);

% [row,colm]=find(OBV==I)



T_Final=T(Andis);
disp('Temperature, T');disp(T_Final)
RH=Delta_Enthalpy(Andis);
RS=Delta_Entropy(Andis);


plot(T,Delta_Enthalpy,T,Delta_Entropy,T,Delta_Internal_Energy,T,Delta_Helmholtz_Energy)

legend({'Residual Enthalpy','Residual Entropy',...
    'Residual Internal Energy',...
    'Residual Helmholtz Energy'})



