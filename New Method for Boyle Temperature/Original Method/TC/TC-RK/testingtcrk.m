clc
clear
close all

format long

P=1e-07;R=83.144598;

Tc=190.564;Pc=45.99;w=0.0115;Zra=0.2894;L=0.217;M=0.9082;N=1.8172; %Methane
% Tc=126.2;Pc=34;w=0.0377;Zra=0.2908;L=0.1901;M=0.89;N=2.0107; %N2
% Tc=150.86;Pc=48.98;w=0;Zra=0.2941;L=0.202;M=0.9085;N=1.8138; %Ar
% Tc=132.92;Pc=34.99;w=0.0482;Zra=0.2903;L=0.1623;M=0.8778;N=2.1566; %CO
% Tc=304.21;Pc=73.83;w=0.2236;Zra=0.2728;L=0.2806;M=0.8684;N=2.2782; %CO2
% Tc=388.37;Pc=27.775;w=0.3558;Zra=0.2719;L=0.2295;M=0.8303;N=2.8936; %RC-318
% Tc=562.05;Pc=48.95;w=0.2103;Zra=0.2698;L=0.1919;M=0.8469;N=2.5997; % benzene
% Tc=405.65;Pc=112.8;w=0.2526;Zra=0.246;L=0.2988;M=0.8652;N=2.3219; % ammonia
% Tc=512.5;Pc=80.84;w=0.5658;Zra=0.2263;L=0.7619;M=0.9142;N=1.7582; % methanol
% Tc=647.096;Pc=220.64;w=0.3449;Zra=0.2289;L=0.4171;M=0.8758;N=2.1818; % water







a=(1/(9*((2^(1/3))-1)))*((((R^2)*(Tc^2)))/Pc);
   
    
            c=((R*Tc)/Pc)*(0.2150-(0.7314*Zra));
        
        b=((((2^(1/3))-1)/3)*((R*Tc)/Pc))-c;

        
        x1=input('Enter lowe range of temperature, T1=');
x2=input('Enter uper range of temperature, T2=');
T=linspace(x1,x2,2000)';
% % % T=     4.422443561780891e+02;


l1=numel(T);
         x0=zeros(l1,1);
for i=1:l1
    
    
    x0(i)=(R*T(i))/P;
    
    
end
        
V=zeros(l1,1);Tr=zeros(l1,1);
alpha=zeros(l1,1);A11=zeros(l1,1);B1=zeros(l1,1);C1=zeros(l1,1);

        for i=1:l1
            Tr(i)=T(i)/Tc;
            
        alpha(i)=(Tr(i)^(N*(M-1)))*(exp(L*(1-(Tr(i)^(M*N)))));
        
        f=@(x)((R*T(i))/(x-b))-((a*alpha(i))/(((x+c)*(x+b+(2*c)))))-P;

        V(i)=fzero(f,x0(i)); % OF : Volume




        end


        
         B=zeros(l1,1);C=zeros(l1,1); OBJ=zeros(l1,1);
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
        
      H_R(i)=(R*T(i)*(Z(i)-1))-((A1(i)/((b+c)))*log((Z(i)+C(i))/...
        (Z(i)+(2*C(i))+B(i))));
    
    A2(i)=a*first_order_alpha(i);
    
%     S_R(i)=(R*(log(Z(i)-B(i))))+((A2(i)/((2^(3/2))*(b+c)))*log((Z(i)+(C(i)*(2+sqrt(2)))+...
%         (B(i)*(1+sqrt(2))))/(Z(i)+(C(i)*(2-sqrt(2)))+(B(i)*(1-sqrt(2))))));
%   
        S_R(i)=(R*(log(Z(i)-B(i))))+((A2(i)/((b+c)))*log(((Z(i)+(2*C(i))+B(i)))/(Z(i)+C(i))));

    F_B(i)=H_R(i)/S_R(i);
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

%    B1=zeros(l1,1);ax1=zeros(l1,1);
%         for i=1:l1        
%     B1(i)=b-((a*alpha(i))/(R*T(i)));
%     ax1(i)=((B1(i)*Pc)/(R*Tc));
%         end
%         
%         
%         KM=[ax1 F_B];



