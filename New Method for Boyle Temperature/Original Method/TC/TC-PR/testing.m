clc
clear
close all

format long

P=1e-07;R=83.144598;
% Tc=510;Pc=60.8;w=0.1986;Zra=0.2633;L=1.3448;M=1;N=0.4642; % R30
% Tc=514;Pc=61.37;w=0.6436;Zra=0.2377;L=0.777;M=0.9043;N=1.6691; % Ethanol
% Zra=0.2941;L=0.1227;M=0.9045;N=1.8541;Tc=150.86;Pc=48.98;w=0; %Argon
% Zra=0.2289;L=0.3872;M=0.8720;N=1.9668;Tc=647.096;Pc=220.64;w=0.3449; %Water
% Tc=487.2;Pc=60;w=0.2556;Zra=0.257;L=0.782;M=0.8656;N=0.8912; 
% Tc=126.2;Pc=34;w=0.0377;Zra=0.2908;L=0.1242;M=0.8898;N=2.013; %N2
% Tc=132.92;Pc=34.99;w=0.0482;Zra=0.2903;L=0.0997;M=0.8782;N=2.1509; %CO
% Tc=611.3;Pc=34.46;w=0.5586;Zra=0.2515;L=1.5803;M=1;N=0.703; %hexanol
% Tc=190.564;Pc=45.99;w=0.0115;Zra=0.2894;L=0.1473;M=0.9075;N=1.8243; %methane
% Tc=304.21;Pc=73.83;w=0.2236;Zra=0.2728;L=0.1784;M=0.859;N=2.4107; % CO2
% Tc=562.05;Pc=48.95;w=0.2103;Zra=0.2698;L=0.1348;M=0.8481;N=2.5791; % benzene
% Tc=405.65;Pc=112.8;w=0.2526;Zra=0.246;L=0.2274;M=0.8645;N=2.3318; % ammonia
Tc=512.5;Pc=80.84;w=0.5658;Zra=0.2263;L=0.665;M=0.9116;N=1.7833; %methanol




a=(0.45724*(R^2)*(Tc^2))/Pc;
c=((R*Tc)/Pc)*(0.1975-(0.7325*Zra));
        b=((0.07780*R*Tc)/Pc)-c;
%         
x1=input('Enter lowe range of temperature, T1=');
x2=input('Enter uper range of temperature, T2=');
T=linspace(x1,x2,1000)';
% T=     4.774880220110055e+02;


l1=numel(T);
         x0=zeros(l1,1);
for i=1:l1
    
    
    x0(i)=(R*T(i))/P;
    
    
end
        
V=zeros(l1,1);Tr=zeros(l1,1);
alpha=zeros(l1,1);A11=zeros(l1,1);B1=zeros(l1,1);C1=zeros(l1,1);Delta_Gibbs_free_Energy=zeros(l1,1);

        for i=1:l1
            Tr(i)=T(i)/Tc;
            
        alpha(i)=(Tr(i)^(N*(M-1)))*(exp(L*(1-(Tr(i)^(M*N)))));
%         
%         f=@(x)((R*T(i))/(x-b))-((a*alpha(i))/(((x+c)*(x+b+(2*c)))+((b+c)*(x-b))))-P;
% 
%         V(i)=fzero(f,x0(i)); % OF : Volume

A11(i)=((a*alpha(i)*P)/((R*T(i))^2));
B1(i)=((b*P)/(R*T(i)));
C1(i)=((c*P)/(R*T(i)));

    X=roots([-1 (-(4*C1(i))-B1(i)+1) ((4*C1(i))+(2*B1(i))-A11(i)-(2*(C1(i)^2))+(3*(B1(i)^2))+(4*B1(i)*C1(i)))...
        ((2*(C1(i)^2))-(B1(i)^2)+(A11(i)*B1(i))-(B1(i)^3)+(2*B1(i)*(C1(i)^2)))]);

    Zv=max(real(X));Zl=min(real(X));
    z=[Zv Zl];

V(i)=(max(z)*R*T(i))/P;



        end
        
  
        B=zeros(l1,1);C=zeros(l1,1);OBJ=zeros(l1,1);
        Z=zeros(l1,1);A1=zeros(l1,1);
        A2=zeros(l1,1);first_order_alpha=zeros(l1,1);
        H_R=zeros(l1,1);S_R=zeros(l1,1);F_B=zeros(l1,1);
        
        for i=1:l1
            
    B(i)=(b*P)/(R*T(i));
     C(i)=(c*P)/(R*T(i));
    Z(i)=(P*V(i))/(R*T(i));
    
    first_order_alpha(i)=(((M-1)*(T(i)^((N*(M-1))-1))*N*(exp(L*(1-(Tr(i)^(N*M))))))/(Tc^(N*(M-1))))-...
    (L*N*(exp(L*(1-(Tr(i)^(N*M)))))*M*(T(i)^((2*N*M)-N-1)))/(Tc^(N*(-1+(2*M))));

       A1(i)=a*((first_order_alpha(i)*T(i))-alpha(i));


    H_R(i)=(R*T(i)*(Z(i)-1))+((A1(i)/((2^(3/2))*(b+c)))*log((Z(i)+(C(i)*(2+sqrt(2)))+...
        (B(i)*(1+sqrt(2))))/(Z(i)+(C(i)*(2-sqrt(2)))+(B(i)*(1-sqrt(2))))));
    
    A2(i)=a*first_order_alpha(i);
    
    S_R(i)=(R*(log(Z(i)-B(i))))+((A2(i)/((2^(3/2))*(b+c)))*log((Z(i)+(C(i)*(2+sqrt(2)))+...
        (B(i)*(1+sqrt(2))))/(Z(i)+(C(i)*(2-sqrt(2)))+(B(i)*(1-sqrt(2))))));
  
    
    F_B(i)=H_R(i)/S_R(i);
    
        Delta_Gibbs_free_Energy(i)=((H_R(i)-(T(i)*S_R(i)))/(S_R(i)));

    OBJ(i)=x0(i)-V(i);
    
        end
        
        
        KM=[T F_B];
        
        I=max(F_B);
disp('maximum of Residual Enthalpy per Residual Entropy=');disp(I);

Andis=find(F_B==I);
% [row,colm]=find(OBV==I)

T_Final=T(Andis);
disp('The Boyle temperature, T');disp(T_Final)


% plot(T_Final,I,'ok',T,F_B,'r')

plot(T,F_B,'r')
% legend({'Boyle Temperature','Ratio of Residual Enthalpy to Residual Entropy'})

        
%         
%         B1=zeros(l1,1);ax1=zeros(l1,1);
%         for i=1:l1        
%     B1(i)=b-((a*alpha(i))/(R*T(i)));
%     ax1(i)=((B1(i)*Pc)/(R*Tc));
%         end
%         
%         
%         KM=[ax1 F_B];
%         
%         
x=zeros(1000,1);
for i=1:1000
    x(i)=T_Final;
end
y=linspace(1,I,1000)';


MP=[T F_B T Delta_Gibbs_free_Energy x y];