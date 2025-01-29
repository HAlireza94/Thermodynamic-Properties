clc
clear
close all

format long

N=input('Enter Group, name=','s');
name=lower(N);
switch name
    case 'alkanols'
        data=load('Alkanols');
    case 'thiophenes'
        data=load('Thiophenes');
    case 'pyridines'
        data=load('Pyridines');
    case 'alkanes'
        data=load('Alkanes');
    case 'alkenes'
        data=load('Alkenes');
    case 'cycloalkanes'
        data=load('Cycloalkanes');
    case 'amines'
        data=load('Amines');
    case 'glycol ethers'
        data=load('Glycol ethers');
    case 'water'
        data=load('Water');
    case 'aromatics'
        data=load('Aromatics');
    case 'gases'
        data=load('Gases');
    case 'ethers'
        data=load('Ethers');
    case 'ketones'
        data=load('Ketones');
         case 'halogenes'
        data=load('Halogenes');
        case 'noble gases'
        data=load('Noble gases');
        case 'polar gases'
        data=load('Polar gases');
end

Tc=data.Tc;
pc=data.Pc;
w=data.w;

n1=numel(Tc);


P=1e-2;R=8.3144598;

Pc=pc.*10^5;

b=zeros(n1,1);
a=zeros(n1,1);
for i=1:n1
    
 a(i)=(((0.42747*(R^2)*(Tc(i)^2)))/Pc(i));
  b(i)=(0.08664*R*Tc(i))/Pc(i);
  
end




T_B=[304.95
488
315.79
616.96
1241.89
862.77];

M=input('Enter Name of Alpha function, M=','s');

switch lower(M)
    case 'nasrifar'
        
        X=zeros(n1,1);
        T_inv=zeros(n1,1);
        k=zeros(n1,1);
        for i=1:n1
 k(i)=0.480+(1.574*w(i))-(0.176*w(i)^2);
  b1(i)=0.25*(12-(11*k(i))+(k(i)^2));
  b2(i)=0.5*(-6+(9*k(i))-(k(i)^2));
  b3(i)=0.25*(4-(7*k(i))+(k(i)^2));
  
            F=@(x)   -(((((a(i)*(-((Tc(i)*b1(i))/(x^2))-...
        ((2*(Tc(i)^2)*b2(i))/(x^3))-...
        ((3*(Tc(i)^3)*b3(i))/(x^4))))*R*x)-...
                (R*a(i)*((((b1(i)/(x/Tc(i)))+...
         (b2(i)/((x/Tc(i)))^2)+(b3(i)/((x/Tc(i)))^3)))))))/...
    ((R*x)^2))-((b(i)-((a(i)*((((b1(i)/(x/Tc(i)))+...
         (b2(i)/((x/Tc(i)))^2)+(b3(i)/((x/Tc(i)))^3)))))/(R*x)))/x);

T_inv(i)=fzero(F,T_B(i));

X(i)=T_inv(i)/T_B(i);

% T_inv_r(i)=T_inv(i)/Tc(i);
% T_B_r(i)=T_B(i)/Tc(i);
% Pr(i)=P/Pc(i);

        end
        
    case 'original'

                X=zeros(n1,1);
        T_inv=zeros(n1,1);
        k=zeros(n1,1);
        for i=1:n1
           k(i)=0.480+(1.574*w(i))-(0.176*w(i)^2); %original 


            F=@(x)  -(((((a(i)*(-(k(i)/x)*(sqrt(x/Tc(i)))*...
                (sqrt(((1+(k(i)*(1-((x/Tc(i))^0.5))))^2)))))*R*x)-...
                (R*a(i)*((1+(k(i)*(1-((x/Tc(i))^0.5))))^2))))/...
    ((R*x)^2))-((b(i)-((a(i)*((1+(k(i)*(1-((x/Tc(i))^0.5))))^2))/(R*x)))/x);

T_inv(i)=fzero(F,T_B(i));

X(i)=T_inv(i)/T_B(i);

% T_inv_r(i)=T_inv(i)/Tc(i);
% T_B_r(i)=T_B(i)/Tc(i);
% Pr(i)=P/Pc(i);

        end

end





