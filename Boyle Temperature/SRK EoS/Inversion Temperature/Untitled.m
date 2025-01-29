clc
clear
close all

%     Tc=126.2;Pc=34;w=0.039;Zc=0.294; % Nitrogrn
%     Tc=304.2;Pc=73.8;w=0.239; % Carbon dioxide
%     Tc=190.6;Pc=46.1;w=0.011; % Methane
        Tc=305.3;Pc=49;w=0.099;Zc=0.279; % Ethane




    R=83.14472;a=(((0.42747*(R^2)*(Tc^2)))/Pc);
b=(0.08664*R*Tc)/Pc;
k=0.480+(1.574*w)-(0.176*w^2);

T_B=693.54;

G=@(x) -(((((a*(-(k/x)*(sqrt(x/Tc))*...
        (sqrt(((1+(k*(1-((x/Tc)^0.5))))^2)))))*R*x)-...
        (R*a*((1+(k*(1-((x/Tc)^0.5))))^2))))/...
    ((R*x)^2))-((b-((a*((1+(k*(1-((x/Tc)^0.5))))^2))/(R*x)))/x);



T_inv=fzero(G,T_B);


T=linspace(0.5*Tc,T_inv-0.05,1000)';
p=Pc;

 Pr=zeros(1000,1);
         Tr=zeros(1000,1);

for i=1:1000
    
     x0(i)=(R*T(i))/p;
    Tr(i)=T(i)/Tc;
    
    alpha(i)=((1+(k*(1-(Tr(i)^0.5))))^2);
    first_order_alpha(i)=-(k/T(i))*(sqrt(Tr(i)))*(sqrt(alpha(i)));
    
    F=@(x) (T(i)*((R/(x-b))-((a*first_order_alpha(i))/(x*(x+b)))))...
        +(x*((((a*alpha(i)*((2*x)+b))/((x*(x+b))^2)))-...
        ((R*T(i))/((x-b)^2))));
    
    
     V(i)=fzero(F,x0(i));
    
    P(i)=((R*T(i))/(V(i)-b))-((a*alpha(i))/(V(i)*(V(i)+b))); %$$$$
    
    Pr(i)=P(i)/Pc;
   
    
    
end
    



