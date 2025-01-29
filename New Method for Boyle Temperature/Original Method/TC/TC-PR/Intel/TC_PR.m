clc
clear
close all
format long

 data=load('SS2');
        
        Tc=data.Tc;Pc=data.Pc;
        w=data.w;Zra=data.Zra;
        L=data.L;M=data.M;N=data.N;
        
        P=1e-07;R=83.144598;

        
        j1=numel(Tc);
        a=zeros(j1,1);b=zeros(j1,1);c=zeros(j1,1);
        for i=1:j1
            
            a(i)=(0.45724*(R^2)*(Tc(i)^2))/Pc(i);
            c(i)=((R*Tc(i))/Pc(i))*(0.1975-(0.7325*Zra(i)));
            b(i)=((0.07780*R*Tc(i))/Pc(i))-c(i);
 
                        
        end
        
        
        x1=input('Enter lowe range of temperature, T1=');
x2=input('Enter uper range of temperature, T2=');
T=linspace(x1,x2,5000)';l1=numel(T);
         x0=zeros(l1,1);
for i=1:l1
        
    x0(i)=(R*T(i))/P;
        
end

Tr=zeros(l1,j1);alpha=zeros(l1,j1);
OF=zeros(l1,j1);

for i=1:j1
    for j=1:l1
        
        Tr(j,i)=T(j)/Tc(i);
        
               alpha(j,i)=(Tr(j,i)^(N(i)*(M(i)-1)))*(exp(L(i)*(1-(Tr(j,i)^(M(i)*N(i))))));
        
        f=@(x)((R*T(j))/(x-b(i)))-((a(i)*alpha(j,i))/(((x+c(i))*(x+b(i)+(2*c(i))))+((b(i)+c(i))*(x-b(i)))))-P;
        
        
          OF(j,i)=fzero(f,x0(j));
   
        
        
    end
    
end


B=zeros(l1,j1);C=zeros(l1,j1);
Z=zeros(l1,j1);first_order_alpha=zeros(l1,j1);
A1=zeros(l1,j1);A2=zeros(l1,j1);H_R=zeros(l1,j1);
S_R=zeros(l1,j1);F_B=zeros(l1,j1);

for i=1:j1
    
    for j=1:l1
        
        
B(j,i)=(b(i)*P)/(R*T(j));
     C(j,i)=(c(i)*P)/(R*T(j));
    Z(j,i)=(P*OF(j,i))/(R*T(j));
    
    first_order_alpha(j,i)=(((M(i)-1)*(T(j)^((N(i)*(M(i)-1))-1))*N(i)*...
        (exp(L(i)*(1-(Tr(j,i)^(N(i)*M(i)))))))/(Tc(i)^(N(i)*(M(i)-1))))-...
    (L(i)*N(i)*(exp(L(i)*(1-(Tr(j,i)^(N(i)*M(i))))))*M(i)*(T(j)^((2*N(i)*M(i))-N(i)-1)))/(Tc(i)^(N(i)*(-1+(2*M(i)))));

       A1(j,i)=a(i)*((first_order_alpha(j,i)*T(j))-alpha(j,i));


    H_R(j,i)=(R*T(j)*(Z(j,i)-1))+((A1(j,i)/((2^(3/2))*(b(i)+c(i))))*log((Z(j,i)+(C(j,i)*(2+sqrt(2)))+...
        (B(j,i)*(1+sqrt(2))))/(Z(j,i)+(C(j,i)*(2-sqrt(2)))+(B(j,i)*(1-sqrt(2))))));
    
    A2(j,i)=a(i)*first_order_alpha(j,i);
    
    S_R(j,i)=(R*(log(Z(j,i)-B(j,i))))+((A2(j,i)/((2^(3/2))*(b(i)+c(i))))*log((Z(j,i)+(C(j,i)*(2+sqrt(2)))+...
        (B(j,i)*(1+sqrt(2))))/(Z(j,i)+(C(j,i)*(2-sqrt(2)))+(B(j,i)*(1-sqrt(2))))));
  
    
    F_B(j,i)=H_R(j,i)/S_R(j,i);
            
                             
        
    end
    
end



I=zeros(j1,1);
for i=1:j1
    
  I(i)=max(F_B(:,i));

end

Andis=zeros(j1,1);
for i=1:j1
    
  Andis(i)=find(F_B(:,i)==I(i));

end

T_Final=zeros(j1,1);
for i=1:j1

    T_Final(i)=T(Andis(i));
    
end