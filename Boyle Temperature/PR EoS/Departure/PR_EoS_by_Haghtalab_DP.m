clc
clear 
close all

format long

% Tc=150.86;Pc=48.98;w=0.000;Zc=0.291; % Argon

% Tc=190.6;Pc=46.1;w=0.011; % Methane

% Tc=305.3;Pc=49;w=0.099; % Ethane

% Tc=369.9;Pc=42.5;w=0.153;  % Propane

% Tc=425.12;Pc=37.96;w=0.2;Zc=0.274; % Butane

% Tc=469.7;Pc=33.7;w=0.252; % Pentane

% Tc=507.6;Pc=30.25;w=0.301;Zc=0.266; % Hexane

% Tc=540.2;Pc=27.4;w=0.350; % Heptane

% Tc=568.7;Pc=24.9;w=0.400;Zc=0.256; % Octane

% Tc=594.6;Pc=22.9;w=0.444;Zc=0.235; % Nonane

% Tc=617.7;Pc=21.1;w=1.5; % Decane w=0.492

% Tc=304.2;Pc=73.8;w=0.239;Zc=0.274; % Carbon dioxid

% Tc=562;Pc=48.9;w=0.212;Zc=0.268; % Benzene

% Tc=513;Pc=81;w=0.556;Zc=0.224; % Methanol

% Tc=647.13;Pc=220.55;w=0.34449;Zc=0.229; %Water

Tc=405.65;Pc=112.8;w=0.2526;Zc=0.242; %Ammonia

%     Tc=1393.24;Pc=1.56;w=2.436;Zc=0.086; % hectane



P=0.0001;
R=83.14472;
a=(0.45724*(R^2)*(Tc^2))/Pc;b=(0.07780*R*Tc)/Pc;

n2=1.7309+(1.6571*w)+(0.1554*w^2); % Haghtalab et al. (2010)

% x1=input('Enter lowe range of temperature, T1=');
% x2=input('Enter uper range of temperature, T2=');


 T=linspace(Tc,1050,100)';n=numel(T);

 
 x0=zeros(n,1);
for i=1:n
    
    
    x0(i)=(R*T(i))/P;
    
    
end


OF=zeros(n,1);
for i=1:n
    
     f=@(x) ((R*T(i))/(x-b))-((a*(exp(1-(n2^(log(T(i)/Tc))))))/((x^2)+(2*b*x)-(b^2)))-P;
     
     OF(i)=fzero(f,x0(i)); % OF : Volume
%           OF(i)=fsolve(f,x0(i));
% OF(i)=fminsearch(f,x0(i))

end
alpha=zeros(n,1);Tr=zeros(n,1);
DEI=zeros(n,1);B=zeros(n,1);
Delta_Gipss=zeros(n,1);A=zeros(n,1);
Delta_Helmholtz_Energy=zeros(n,1);
Delta_Internal_Energy=zeros(n,1);
Delta_Entropy=zeros(n,1);
Delta_Enthalpy=zeros(n,1);
Z=zeros(n,1);
  F_B=zeros(n,1);

for i=1:n
    
    Tr(i)=T(i)/Tc;
    B(i)=(b*P)/(R*T(i));
    Z(i)=(P*OF(i))/(R*T(i));
    
    alpha(i)=(exp(1-(n2^(log(T(i)/Tc))))); % Haghtalab et al. (2010)
    A(i)=(a*alpha(i))/((R*T(i))^2);
    
    
    
    DEI(i)=(a*alpha(i))-((-T(i)*a*alpha(i)*(log(n2))*(Tr(i)^(log(n2))))/T(i));
    

    
    Delta_Enthalpy(i)=((((DEI(i))/(2*sqrt(2)*b)))*...
        (log((OF(i)+((1-sqrt(2))*b))/(OF(i)+((1+sqrt(2))*b)))))+(R*T(i)*(Z(i)-1)); 
    
    Delta_Entropy(i)=(R*(log(Z(i)-B(i))))-...
        (((((-a*alpha(i)*(log(n2))*(Tr(i)^(log(n2))))/T(i)))/(2*sqrt(2)*b))*...
        (log((OF(i)+((1-sqrt(2))*b))/(OF(i)+((1+sqrt(2))*b))))); 
    
    
    Delta_Internal_Energy(i)= Delta_Enthalpy(i)-(R*T(i)*(Z(i)-1));
    
%     Delta_Helmholtz_Energy(i)=((R*T(i))*(log(1/(Z(i)-B(i)))))-(((a*alpha(i))/(2*sqrt(2)*b))*...
%         (log((OF(i)+((1-sqrt(2))*b))/(OF(i)+((1+sqrt(2))*b)))));
    
    
    
    F_B(i)=Delta_Enthalpy(i)/Delta_Entropy(i);
    
    Delta_Gipss(i)=Delta_Enthalpy(i)-(T(i)*Delta_Entropy(i));
    
    Delta_Helmholtz_Energy(i)=Delta_Internal_Energy(i)-(T(i)*Delta_Entropy(i));
end




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
