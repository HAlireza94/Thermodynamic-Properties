clc
clear 
close all



format long

% N=input('Enter Group, name=','s');
% name=lower(N);
% switch name
%     case 'alkanols'
%         data1=load('Alkanols');
%     case 'thiophenes'
%         data1=load('Thiophenes');
%     case 'pyridines'
%         data1=load('Pyridines');
%     case 'alkanes'
%         data1=load('Alkanes');
%         data2=load('Coefficient_Cp');
%         case 'polar gases'
%         data1=load('Polar gases');
%     case 'alkenes'
%         data1=load('Alkenes');
%                 data2=load('Coefficient_Cp_Alkenes');
% 
%     case 'cycloalkanes'
%         data1=load('Cycloalkanes');
%     case 'amines'
%         data1=load('Amines');
%     case 'glycol ethers'
%         data1=load('Glycol ethers');
%     case 'water'
%         data1=load('Water');
%     case 'aromatics'
%         data1=load('Aromatics');
%     case 'gases'
%         data1=load('Gases');
%     case 'ethers'
%         data1=load('Ethers');
%     case 'ketones'
%         data1=load('Ketones');
%         case 'halogens'
%         data1=load('Halogenes');
%         case 'noble gases'
%         data1=load('Noble gases');
%         
%     case 'fatty acid esters'
%         
%         data1=load('Fatty acid esters');
% end

    Tc=190.56;pc=45.99200;w=0.011;Zc=0.2863;


% Tc=data1.Tc;
% pc=data1.Pc;
% w=data1.w;
% Zc=data1.Zc;

% M=data1.Mw;

M=16.043;

Pc=pc.*10^5;

n1=numel(Tc);

P=1e-2;R=8.3144598;

% A=data2.A;
% B=data2.B;
% C=data2.C;
% D=data2.D;
% E=data2.E;

A=34.942;B=-0.039957;C=0.00019184;D=-1.5303e-07;E=3.9321e-11;

for i=1:n1
    
    a(i)=(((0.42747*(R^2)*(Tc(i)^2)))/Pc(i));
  b(i)=(0.08664*R*Tc(i))/Pc(i);
  
  c1(i)=(-45.7247*((1/3)-Zc(i)));
  c2(i)=((-2.184*exp(c1(i)))+0.2658);
  c(i)=((((1/3)-Zc(i))*((R*Tc(i))/Pc(i)))*c2(i));
  
  m(i)=0.266+(0.4459*w(i)^0.5);
  n4(i)=(1/m(i))*(0.2469+(0.7495*w(i)));

end

Pr=P/Pc;



T=[90.694
100.69
110.69
120.69
130.69
140.69
150.69
160.69
170.69
180.69
190.69
200.69
210.69
220.69
230.69
240.69
250.69
260.69
270.69
280.69
290.69
300.69
310.69
320.69
330.69
340.69
350.69
360.69
370.69
380.69
390.69
400.69
410.69
420.69
430.69
440.69
450.69
460.69
470.69
480.69
490.69
500.69
510.69
520.69
530.69
540.69
550.69
560.69
570.69
580.69
590.69
600.69
610.69
620.69];


 n=numel(T);

 
 x0=zeros(n,1);
for i=1:n
    
    
    x0(i)=(R*T(i))/P;
    
    
end


OF=zeros(n,1);
for i=1:n
    
     f=@(x) ((R*T(i))/(x+c-b))-((a*((exp(m*(1-((T(i)/Tc)^n4))))^2))/((x+c)*(x+c+b)))-P;
     
     OF(i)=fzero(f,x0(i)); % OF : Volume
%           OF(i)=fsolve(f,x0(i));
% OF(i)=fminsearch(f,x0(i))

end


J_T=zeros(n,1);
Tr=zeros(n,1);
Cv=zeros(n,1);
Cp=zeros(n,1);
Gamma=zeros(n,1);
X=zeros(n,1);
RSS=zeros(n,1);
sound_speed=zeros(n,1);

for i=1:n
    
    Cp_star(i)=A+(B*T(i))+(C*(T(i)^2))+(D*(T(i)^3))+(E*(T(i)^4));
    
    Cv_star(i)=Cp_star(i)-R;
    
    Tr(i)=T(i)/Tc;
    
    Z(i)=(P*OF(i))/(R*T(i));
    
    alpha(i)=((exp(m*(1-((T(i)/Tc)^n4))))^2);
    
    first_order_alpha(i)=((-2*m*n4)/Tc)*(Tr(i)^(n4-1))*((exp(2*m*(1-((T(i)/Tc)^n4)))));
     

    second_order_alpha(i)=(((-2*m*n4)*(T(i)^(n4-2))*...
        ((exp(2*m*(1-((T(i)/Tc)^n4)))))*...
        (((Tc^n4)*(n4-1))-(2*m*n4*(T(i)^n4))))/...
        (Tc^(2*n4)));
   
     a_first_order_alpha(i)=a*first_order_alpha(i);

     a_second_order_alpha(i)=a*second_order_alpha(i);


    Cv(i)=Cv_star(i)+(T(i)*...
     ((a_second_order_alpha(i))/(b))*...
     (log((b/(OF(i)+c))+1)));
    
    first_order_P_T(i)=(R/(OF(i)+c-b))-...
        (a_first_order_alpha(i)/((OF(i)+c)*(OF(i)+c+b)));
    
    first_order_P_V(i)=((-R*T(i))/((OF(i)+c-b)^2))+...
        ((a*alpha(i)*((2*c)+(2*OF(i))+b)))/((((OF(i)+c)*(OF(i)+c+b)))^2);
    
    Cp(i)=Cv(i)-((T(i)*(first_order_P_T(i)^2))/(first_order_P_V(i)));
    
    X(i)=Cp(i)*Cv(i);
    
    Gamma(i)=Cp(i)/Cv(i);
  
    SS(i)=sqrt((-OF(i)^2)*(Gamma(i)/M)*first_order_P_V(i)); 
    
    sound_speed(i)=SS(i)*sqrt(1000); %$#@
    
    Kappa(i)=((-(1/OF(i)))*(first_order_P_V(i)^(-1)));
    
    
    J_T(i)=((T(i)*((-first_order_P_T(i))/first_order_P_V(i)))-OF(i))/Cp(i);
      
    
end
    
    

figure (1)

plot(T,Cp)
xlabel('Temperature (K)');ylabel('heat capacity at constant pressure(Cp)')

figure (2)

plot(T,Cv)
xlabel('Temperature (K)');ylabel('heat capacity at constant volume (Cv)')


figure (3)
plot(T,Gamma)
xlabel('Temperature (K)');ylabel('Gamma (Cp/Cv)')




