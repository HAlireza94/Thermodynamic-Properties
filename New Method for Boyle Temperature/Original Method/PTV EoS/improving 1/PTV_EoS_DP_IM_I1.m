clc
clear 
close all

format long

        data=load('SS1');
        
Tc=data.Tc;Pc=data.Pc;
w=data.w;Zc=data.Zc;

n1=numel(Tc);
P=1e-7;R=83.14472;
F=zeros(n1,1);b=zeros(n1,1);a=zeros(n1,1);c=zeros(n1,1);
Oa=zeros(n1,1);Ob=zeros(n1,1);Oc=zeros(n1,1);
for i=1:n1
    

   Oa(i)=(0.66121-(0.76105*Zc(i)));
   Ob(i)=(0.02207+(0.20868*Zc(i)));
   Oc(i)=(0.57765-(1.87080*Zc(i)));

    a(i)=(((Oa(i)*(R^2)*(Tc(i)^2)))/Pc(i));
   b(i)=((Ob(i)*R*Tc(i))/Pc(i));
 c(i)=(Oc(i)*R*Tc(i))/Pc(i);
    F(i)=0.46283+(3.58230*w(i)*Zc(i))+(8.19417*(w(i)^2)*(Zc(i)^2));
    
end


x1=input('Enter lowe range of temperature, T1=');
x2=input('Enter uper range of temperature, T2=');


 T=linspace(x1,x2,2500)';n=numel(T);
 
 
 x0=zeros(n,1);
for i=1:n
    
    
    x0(i)=(R*T(i))/P;
    
    
end

T_B_Tsonopoulos=[503.26 770.15 908.1 1019.23 1097.94 1159.7 ...
1205.78 1241.45 1271.78 1295.73 989.47 1108.69 1229.83 ...
1322.53 1340.91 1375.28 1407.34 717 899.86 1012.45 1094.92 1159.68 ...
1202.63 1241.64 1138.37 1170.39 1219.32 1141.13 1159.99 ...
1158.45 1215.45 1255.52 963.36 1031.54 1107.13 328.75 ...
400.75 344.62 720.25 947.95 1339.67 551.29 ...
1268.75 1223.34 812.25 1022.71 761.81 864.8 ...
1037.94 1006.56 1060.83 1125.32 1185.54 1228.66 ...
992.7 1042.71 1149.48 1194.68 1237.1 1449.11];

% T_B_Meng_Duan=[510.39 757.52 880.03 ...
% 976.45 1040.84 1090.45 1125.77 1152.50 1175.57 1193.65 965.33 1064.79 ...
% 1179.53 1264.81 1281.85 1300.92 1322.48 707.46 874.97 967.69 1041.14 ...
% 1095.85 1123.70 1153.53 1078.84 1103.16 1143.33 1068.05 1080.77 ...
% 1086.77 1137.60 1172.09 924.27 980.66 1047.11 330.12 408.22 ...
% 345.06 686.38 898.54 1284.18 530.54 1208.34 1172.38 729.36 ...
% 989.86 731.28 806.52 898.77 884.11 944.13 1011.41 1070.82 ...
% 1115.30 946.11 986.59 1088.22 1122.73 1158.78 1325.29];


OF=zeros(n,n1);
for j=1:n1
    
for i=1:n
    
     f=@(x) ((R*T(i))/(x-b(j)))-((a(j)*((1+(F(j)*(1-((T(i)/Tc(j))^0.5))))^2))/((x*(x+b(j)))+(c(j)*(x-b(j)))))-P;
     
     OF(i,j)=fzero(f,x0(i));
%           OF(i)=fsolve(f,x0(i));
% OF(i)=fminsearch(f,x0(i))

end

end

DEI=zeros(n,1);C=zeros(n,1);
A=zeros(n,1);alpha=zeros(n,1);B=zeros(n,1);Tr=zeros(n,1);
Delta_Gipss=zeros(n,1);Delta_Helmholtz_Energy=zeros(n,1);
Delta_Internal_Energy=zeros(n,1);Delta_Entropy=zeros(n,1);
Delta_Enthalpy=zeros(n,1);Z=zeros(n,1);F_B=zeros(n,1);X21=zeros(n,1);
TSV=zeros(n,1);first_D_alpha=zeros(n,1);second_D_alpha=zeros(n,1);
theird_D_alpha=zeros(n,1);SV=zeros(n,1);X18=zeros(n,1);X20=zeros(n,1);
Error=inf;
betta=linspace(-100,100,1000)';
n2=numel(betta);

for j=1:n1

for i=1:n
    
    Tr(i,j)=T(i)/Tc(j);
    B(i,j)=(b(j)*P)/(R*T(i));
    C(i,j)=(c(j)*P)/(R*T(i));
    Z(i,j)=(P*OF(i,j))/(R*T(i));
    alpha(i,j)=((1+(F(j)*(1-((T(i)/Tc(j))^0.5))))^2);
    
    A(i,j)=(a(j)*alpha(i,j))/((R*T(i))^2);
     
    DEI(i,j)=(a(j)*alpha(i,j)-(T(i)*(((-F(j)*a(j))/T(i))*sqrt(Tr(i,j))*sqrt(alpha(i,j)))));
  
    Delta_Enthalpy(i,j)=(R*T(i)*(Z(i,j)-1))+...
        ((DEI(i,j)/(sqrt((c(j)^2)+(6*b(j)*c(j))+(b(j)^2))))*(log((((2*OF(i,j))-...
        (sqrt((c(j)^2)+(6*b(j)*c(j))+(b(j)^2)))+c(j)+b(j))/((2*OF(i,j))+...
        (sqrt((c(j)^2)+(6*b(j)*c(j))+(b(j)^2)))+c(j)+b(j))))));
    
    Delta_Entropy(i,j)=(R*(log(Z(i,j)-B(i,j))))-...
        (((((-F(j)*a(j))/T(i))*sqrt(Tr(i,j))*sqrt(alpha(i,j))/...
        (sqrt((c(j)^2)+(6*b(j)*c(j))+(b(j)^2)))))*(log((((2*OF(i,j))-...
        (sqrt((c(j)^2)+(6*b(j)*c(j))+(b(j)^2)))+c(j)+b(j))/((2*OF(i,j))+...
        (sqrt((c(j)^2)+(6*b(j)*c(j))+(b(j)^2)))+c(j)+b(j))))));
    
        Delta_Internal_Energy(i,j)=Delta_Enthalpy(i,j)-(R*T(i)*(Z(i,j)-1));
    
    Delta_Helmholtz_Energy(i,j)=((R*T(i))*(log(1/(Z(i,j)-B(i,j)))))-...
        (((a(j)*alpha(i,j))/(sqrt((c(j)^2)+(6*b(j)*c(j))+(b(j)^2))))*(log((((2*OF(i,j))-...
        (sqrt((c(j)^2)+(6*b(j)*c(j))+(b(j)^2)))+c(j)+b(j))/((2*OF(i,j))+...
        (sqrt((c(j)^2)+(6*b(j)*c(j))+(b(j)^2)))+c(j)+b(j))))));
    
      SV(i,j)=b(j)-((a(j)*alpha(i,j))/(R*T(i)));
    
      first_D_alpha(i,j)=(-F(j)/T(i))*sqrt(Tr(i,j))*sqrt(alpha(i,j));
    
      second_D_alpha(i,j)=((F(j)*((F(j)*sqrt(Tr(i,j)))+sqrt(alpha(i,j))))/(2*sqrt(Tc(j))*(T(i)^1.5)));
    
      theird_D_alpha(i,j)=(-0.75*F(j)*sqrt(Tr(i,j))*((F(j)*sqrt(Tr(i,j)))+sqrt(alpha(i,j))))/(T(i)^3);
    
      TSV(i,j)=((((-(T(i)^3)*theird_D_alpha(i,j))+(3*(T(i)^2)*second_D_alpha(i,j))-...
        (6*T(i)*first_D_alpha(i,j))+(6*alpha(i,j)))/(R*T(i)^4)))*a(j);
    
      X18(i,j)=TSV(i,j)*Pc(j)*Tc(j)^2;
                   
    X20(i,j)=-31.331331331331327*X18(i,j);
    
    X21(i,j)=(SV(i,j)*Pc(j))/(X20(i,j));
                
    F_B(i,j)=(Delta_Enthalpy(i,j)/Delta_Entropy(i,j))-X21(i,j);

end


end


% for m=1:n2
%         
%         for j=1:n1
%             
%             for i=1:n
%     
%     
%      X20(i,j,m)=betta(m)*X18(i,j);
%     
%     X21(i,j,m)=(SV(i,j)*Pc(j))/(X20(i,j,m));
%                 
%     F_B(i,j,m)=(Delta_Enthalpy(i,j)/Delta_Entropy(i,j))-X21(i,j,m);
%     
%     
%             end
%         end
%         
%         I=zeros(n1,n2);
%    for i=1:n1
%     
%   I(i,m)=max(F_B(:,i,m));
%   
%    end
%    
%    Andis=zeros(n1,n2);
%    
% for i=1:n1
%     
%   Andis(i,m)=find(F_B(:,i,m)==I(i,m));
% 
% end
% 
% T_Final=zeros(n1,n2);
% 
% for i=1:n1
% 
%     T_Final(i,m)=T(Andis(i,m));
%     
% end
% 
% 
% E=zeros(n1,n2);
% 
% for i=1:n1
%     
%     E(i,m)=(abs(T_Final(i,m)-T_B_Tsonopoulos(i))/T_B_Tsonopoulos(i))*100;
%         
% end
% 
% % Error(m)=(sum(E(:,m)))/60;
% 
% Error(m)=mean(E(:,m));
% 
% end



I=zeros(n1,1);
for i=1:n1
    
  I(i)=max(F_B(:,i));

end

Andis=zeros(n1,1);
for i=1:n1
    
  Andis(i)=find(F_B(:,i)==I(i));

end

T_Final=zeros(n1,1);
for i=1:n1

    T_Final(i)=T(Andis(i));
    
end
