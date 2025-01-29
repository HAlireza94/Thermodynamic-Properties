clc
clear 
close all

Yi=[34.1 33.5 32.4];
Kij=[0 0 0;0 0 0;0 0 0];

Tc=[304.19 190.56 126.2];
Pc=[73.82 45.99 34.6];
w=[0.228 0.012 0.0377];
Zc=[0.274 0.286 0.294];

N=numel(Yi);
pcm=zeros(1,N);
tcm=zeros(1,N);
for i=1:N
    
    pcm(i)=Pc(i)*Yi(i);
    tcm(i)=Tc(i)*Yi(i);

end
PcM=sum(pcm);TcM=sum(tcm);


N=numel(Yi);

    alpha_function=input('Enter name of alpha function :','s');
    
    T1=linspace(150,3000,1000)';
    OBJ=zeros(1000,1);
    for i=1:1000
        
        OBJ(i)=Calculation_T_Inv_VTPR_EoS_for_Mixture(alpha_function,T1(i),Tc,Pc,w,Kij,Yi,N,Zc);

        
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
q=input('Enter the number of devides, q=');
    
    
    T=linspace(150,T_inv,q)';
    
P=zeros(q,1);P1=zeros(q,1);
Pr=zeros(q,1);Tr=zeros(q,1);

for i=1:q
    
P1(i)=Calculation_JTIC_VTPR_EoS_For_Mixture(alpha_function,T(i),Tc,Pc,w,N,Yi,Kij,Zc);
P(i)=0.1*P1(i);
Pr(i)=P1(i)/PcM;
Tr(i)=T(i)/TcM;
end


figure (2)

plot(T,P)
M=[Tr Pr];


% M=[];
% T=M(:,1);px=M(:,2);pex=10*px;
% q=numel(T);pexr=zeros(q,1);
% P=zeros(q,4);Tr=zeros(q,1);Pr=zeros(q,4);
% for i=1:4
%     
%      if i==1
%     
%         alpha_function='original';
%     
%     elseif i==2
%         
%         alpha_function='coquelet';
%         
%     elseif i==3
%         
%         alpha_function='haghtalab';
%         
%     elseif i==4
%         
%         alpha_function='saffari';
% 
%     end
%     
%  for j=1:q
%  
% P(j,i)=Calculation_JTIC_VTPR_EoS_For_Mixture(alpha_function,T(j),Tc,Pc,w,N,Yi,Kij,Zc);
%         
% % P(j,i)=0.1*P1(j,i);
% Pr(j,i)=P(j,i)/PcM;
% Tr(j)=T(j)/TcM;
% pexr(j)=pex(j)/PcM;
%  end
% end
% 
% plot(Tr,Pr(:,1),'k',Tr,Pr(:,2),'r',Tr,Pr(:,3),'b',Tr,Pr(:,4),'g')
% legend({'Original','Coquelet','Haghtalab','Saffari'});
    

