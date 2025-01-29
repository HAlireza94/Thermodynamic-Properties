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
cp=[TcM PcM];
    alpha_function=input('Enter name of alpha function :','s');
    
    
T1=linspace(100,2500,1000)';

 
    OBJ=zeros(1000,1);
    for i=1:1000
        
        OBJ(i)=Calculation_T_Inv_VTSRK_EoS_for_Mixture...
    (alpha_function,T1(i),Tc,Pc,w,Kij,Yi,N,Zc);

        
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

q=input('Enter the number of Devides, q=');

T=linspace(100,T_inv,q)';

P=zeros(q,1);P1=zeros(q,1);
Pr=zeros(q,1);Tr=zeros(q,1);

for i=1:q
    
P1(i)=Calculation_JTIC_VTSRK_EoS_For_Mixture(alpha_function,T(i),Tc,Pc,w,N,Yi,Kij,Zc);
P(i)=0.1*P1(i);
Pr(i)=P1(i)/PcM;
Tr(i)=T(i)/TcM;

end

figure(2)
plot(T,P)

M=[Tr Pr];

% M=[];T=M(:,1);px=M(:,2);
% pex=10*px;
%  q=numel(T); 
% Tr=zeros(q,1); 
% pexr=zeros(q,1); 
% Pr=zeros(q,4);P=zeros(q,4);
% 
% for i=1:4
%     
%      if i==1
%     
%         alpha_Function='original';
%     
%     elseif i==2
%         
%         alpha_Function='ozokwelu';
%         
%     elseif i==3
%         
%         alpha_Function='trebble';
%         
%     elseif i==4
%         
%         alpha_Function='nasrifar';
% 
%     end
%     
%  for j=1:q
%  
% P(j,i)=Calculation_JTIC_VTSRK_EoS_For_Mixture(alpha_Function,T(j),Tc,Pc,w,N,Yi,Kij,Zc);
% 
%  Pr(j,i)=P(j,i)/PcM;
% Tr(j)=T(j)/TcM;
% pexr(j)=pex(j)/PcM;
%  end
% end
% 
% plot(Tr,Pr(:,1),'k',Tr,Pr(:,2),'r',Tr,Pr(:,3),'b',Tr,Pr(:,4),'g')
% legend({'Original','ozokwelu','trebble','nasrifar'});
