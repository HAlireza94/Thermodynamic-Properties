clc
clear
close all

    Tc=190.6;pc=46.1;w=0.011;

% A=34.942;B=-0.039957;C=0.00019184;D=-1.5303e-07;E=3.9321e-11;

% P=1e-2;R=8.3144598;

% Pc=pc.*10^5;


a=(0.45724*(R^2)*(Tc^2))/Pc;b=(0.07780*R*Tc)/Pc;

k=0.41287+(1.34494*w)+(0.00421*w^2);

T=525.53;
V=(R*2*T)/P;

PARA0=[2*T,V];
  

% [X,fval,exitflag,output]=Inve(k,a,b,Tc,R,A,B,C,D,E,PARA0);

F=@(x) -(((((a*(-(k/Tc)*(exp(k*(1-(x/Tc))))))*R*x)-(R*a*(exp(k*(1-(x/Tc)))))))/...
    ((R*x)^2))-((b-((a*exp(k*(1-(x/Tc))))/(R*x)))/x);

T_inv=fzero(F,5*T)




    















