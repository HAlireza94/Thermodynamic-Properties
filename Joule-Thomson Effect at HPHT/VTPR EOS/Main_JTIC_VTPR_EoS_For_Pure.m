clc
clear
close all

format long

N=input('Enter Name Of Gas :','s');

[Tc, Pc, w, Zc]=Critical_Properties(N);


disp('==============================================')

alpha_Function=input('name of Alpha Function, alpha :','s');

disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx')
 
T_B=input('Enter Boyle Temperature, TB=');

disp('<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>')
p=3.5*Pc;

q=input('Enter number of divides, q=');
%for non polar we have to decrease initial pressure 
%but for polar substance we should increase initial pressure 
 
         T_inv=Calculation_T_Inv_VTPR_EoS_For_Pure(alpha_Function,w,Tc,Pc,T_B,Zc);
         T=linspace(0.7*Tc,T_inv,q)';

         Pr=zeros(q,1);
         Tr=zeros(q,1);
         
         for i=1:q
             
         [Pr(i), Tr(i)]=Calculation_JTIC_VTPR_EoS_For_Pure(alpha_Function,p,Pc,Tc,T(i),w,Zc);
%          P(i)=Pc*Pr(i)*1e-05;


         end
         
         plot(Tr,Pr);xlabel('Tr');ylabel('Pr');
         Tr(q)=T_inv/Tc;Pr(q)=1e-07;

         M=[Tr Pr];

% tr=[0.87108
% 0.993031
% 1.11498
% 1.23693
% 1.37631
% 1.48084
% 1.60279
% 1.72474];
% 
%  q=numel(tr); 
% 
%  T=Tc.*tr;
% Pr=zeros(q,4);
% for i=1:4
%     
%      if i==1
%     
%         alpha_Function='original';
%     
%     elseif i==2
%         
%         alpha_Function='coquelet';
%         
%     elseif i==3
%         
%         alpha_Function='haghtalab';
%         
%     elseif i==4
%         
%         alpha_Function='saffari';
% 
%     end
%     
%  for j=1:q
%  
% [Pr(j,i), Tr(j,i)]=Calculation_JTIC_VTPR_EoS_For_Pure(alpha_Function,p,Pc,Tc,T(j),w,Zc);
% 
% 
%  end
% end
% 
% plot(tr,Pr(:,1),'k',tr,Pr(:,2),'r',tr,Pr(:,3),'b',tr,Pr(:,4),'g')
% legend({'Original','Coquelet','Haghtalab','Saffari'});title(N)
%     
%  
%  
% 
% 
