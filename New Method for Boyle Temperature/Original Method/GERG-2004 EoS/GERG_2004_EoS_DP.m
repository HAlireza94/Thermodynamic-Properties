clc
clear
close all

Tc=190.56;Pc=45.99200;w=0.011;Zc=0.2863;RhoC=0.1627; %Methane

% Tc=150.86;Pc=48.98;w=0.000;Zc=0.291;RhoC=0.5356; %Argon

% Tc=126.2;Pc=34.8;w=0.039;Zc=0.294;RhoC=0.3140; % Nitrogen

% Tc=304.2;Pc=73.8;w=0.239;Zc=0.274;RhoC=0.4682; % Carbon dioxid


% it has to be mentioned, P:[bar] T:[K] Rhoc:[g/ml]

Kpol=6;Kexp=18;

R=83.14472;
P=1e-07;

tk=[0.125
1.125
0.375
1.125
0.625
1.5
0.625
2.625
2.75
2.125
2
1.75
4.5
4.75
5
4
4.5
7.5
14
11.5
26
28
30
16];

dk=[1
1
2
2
4
4
1
1
1
2
3
6
2
3
3
4
4
2
3
4
5
6
6
7];

ck=[0
0
0
0
0
0
1
1
1
1
1
1
2
2
2
2
2
3
3
3
6
6
6
6];

nk=[0.573357042
-1.676068752
0.234052918
-0.219473763
0.016369201
0.015004406
0.098990489
0.583827709
-0.747868676
0.300333029
0.209855438
-0.018590151
-0.157825583
0.127167352
-0.032019744
-0.068049729
0.024291413
0.005144045
-0.01908495
0.005522968
-0.004419739
0.040061417
-0.033752086
-0.002512766];



% delta = Rho/RhoC =====> for finding root delta is showed x and Ideal
% delta is showed x0
% To= Tc/T


T=linspace(90,623,1000)';
To=zeros(1000,1);
x0=zeros(1000,1);
for i=1:1000
    
    x0(i)=P/(R*T(i)*RhoC);
    To(i)=Tc/T(i);
end

delta=zeros(1000,1);

for i=1:1000
    
    f=@(x) (x.*R*T(i)*RhoC.*(1+(x.*(sum(nk(1:Kpol).*dk(1:Kpol).*(x.^(dk(1:Kpol)-1)).*(To(i).^(tk(1:Kpol))))+...
        sum(nk((1+Kpol):(Kpol+Kexp)).*(x.^(dk((1+Kpol):(Kpol+Kexp))-1)).*(dk((1+Kpol):(Kpol+Kexp))-...
        (ck((1+Kpol):(Kpol+Kexp)).*(x.^(ck((1+Kpol):(Kpol+Kexp)))))).*...
        (To(i).^(tk((1+Kpol):(Kpol+Kexp)))).*exp(-x.^(ck((1+Kpol):(Kpol+Kexp)))))))))-P;
    
    F=@(X) arrayfun(f,X);
    
    delta(i)=fzero(F,x0(i));
    
    Rho(i)=delta(i)*RhoC;
    
    V(i)=1/Rho(i);
    
    Vi(i)=(R*T(i))/P;
end
    

a1=zeros(Kpol,1000);
ad1=zeros(Kpol,1000);
at1=zeros(Kpol,1000);

for i=1:1000
    
    for j=1:Kpol

     a1(j,i)=(nk(j)*(delta(i)^dk(j))*(To(i)^tk(j)));
    
     ad1(j,i)=(nk(j)*dk(j)*(delta(i)^(dk(j)-1))*(To(i)^(tk(j))));

     at1(j,i)=(nk(j)*tk(j)*(delta(i)^dk(j))*(To(i)^(tk(j)-1)));
    
    end
end


a2=zeros(Kexp,1000);
ad2=zeros(Kexp,1000);
at2=zeros(Kexp,1000);

for i=1:1000
    
    for j=1:Kexp

       a2(j,i)=(nk(j+Kpol)*(delta(i)^dk(j+Kpol))*(To(i)^tk(j+Kpol))*(exp(-delta(i)^ck(j+Kpol))));
    
    ad2(j,i)=(nk(j+Kpol)*(delta(i)^(dk(j+Kpol)-1))*(dk(j+Kpol)-(ck(j+Kpol)*(delta(i)^ck(j+Kpol))))*...
        (To(i)^(tk(j+Kpol)))*(exp(-delta(i)^(ck(j+Kpol)))));
    
    at2(j,i)=(nk(j+Kpol)*tk(j+Kpol)*(delta(i)^dk(j+Kpol))*(To(i)^(tk(j+Kpol)-1))*(exp(-delta(i)^(ck(j+Kpol)))));
    
    end
end

Ad=sum(ad1)+sum(ad2);
At=sum(at1)+sum(at2);
A=sum(a1)+sum(a2);

residual_Enthalpy=zeros(1000,1);
residual_Entropy=zeros(1000,1);
F_B=zeros(1000,1);RV=zeros(1000,1);

for i=1:1000
    
    residual_Enthalpy(i)=R*T(i)*((To(i)*At(i))+(delta(i)*Ad(i)));
    
    residual_Entropy(i)=((To(i)*At(i))-A(i)+(log(1+(delta(i)*Ad(i)))))*R;
    
    F_B(i)=residual_Enthalpy(i)/residual_Entropy(i);
    
    RV(i)=Vi(i)-V(i);
    
end
    
  
  I=max(F_B);
disp('maximum of Residual Enthalpy per Residual Entropy=');disp(I);

Andis=find(F_B==I);
% [row,colm]=find(OBV==I)

T_Final=T(Andis);
disp('The Boyle temperature, T');disp(T_Final)


plot(T_Final,I,'ok',T,F_B)

    
    




