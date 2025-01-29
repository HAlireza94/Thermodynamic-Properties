clc
clear 
close all

yi=[74.2 10 5 3 1.2 1 5.6];Yi=0.01*yi;
Kij=[0 0 0 0 0 0 0
    0 0 0 0 0 0 0
    0 0 0 0 0 0 0
    0 0 0 0 0 0 0
    0 0 0 0 0 0 0
    0 0 0 0 0 0 0
    0 0 0 0 0 0 0];
Cri_P=[190.56	45.99	0.012	0.286
305.32	48.72	0.1	0.279
369.83	42.48	0.152	0.276
425.12	37.96	0.2	0.274
469.7	33.7	0.252	0.27
507.6	30.25	0.301	0.266
540.2	27.4	0.35	0.261];

Tc=Cri_P(:,1);Pc=Cri_P(:,2);
w=Cri_P(:,3);Zc=Cri_P(:,4);

N=numel(Yi);
pcm=zeros(1,N);
tcm=zeros(1,N);
for i=1:N
    
    pcm(i)=Pc(i)*Yi(i);
    tcm(i)=Tc(i)*Yi(i);

end
PcM=sum(pcm);TcM=sum(tcm);
CP=[TcM PcM];

N=numel(Yi);

  T1=linspace(100,2000,3000)';

 
    OBJ=zeros(3000,1);
    for i=1:3000
        
         OBJ(i)=Calculation_T_Inv_PTV_EoS_for_Mixture(T1(i),Tc,Pc,w,Kij,Yi,N,Zc);

        
    end
    
    
    
    figure (1)
    plot(T1,OBJ);xlabel('T');ylabel('OBJ');
grid on
I=min(abs(OBJ));
disp('Objective Value, Inversion Temperature=');disp(I);


Andis=find(abs(OBJ)==I);
% [row,colm]=find(OBV==I)

T_inv=T1(Andis);
disp('Temperature, T');disp(T_inv)  
    
    q=input('Enter devides, q=');

T=linspace(100,T_inv,q)';

% H=[];T=H(:,1);P=H(:,2);P1=10*P;
% q=numel(T);
P1=zeros(q,1);%P=zeros(q,1);
Pr=zeros(q,1);Tr=zeros(q,1);
 
for i=1:q
    
P1(i)=Calculation_JTIC_PTV_EoS_For_Mixture(T(i),Tc,Pc,w,N,Yi,Kij,Zc);
% P(i)=P1(i)*0.1;
Pr(i)=P1(i)/PcM;
Tr(i)=T(i)/TcM;
end

figure (2)
plot(Tr,Pr)
M=[Tr Pr];