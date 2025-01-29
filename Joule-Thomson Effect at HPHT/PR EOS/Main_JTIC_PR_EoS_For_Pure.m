clc
clear
close all

format long

M=input('Enter Name Of Gas :','s');

[Tc, Pc, w]=Critical_Properties(M);


disp('==============================================')

alpha_Function=input('name of Alpha Function, alpha :','s');

disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx')
 
T_B=input('Enter Boyle Temperature, TB=');
q=input('Enter the number of array, q=');

disp('<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>')
p=1.05*Pc;


         T_inv=Calculation_T_Inv_PR_EoS_For_Pure(alpha_Function,w,Tc,Pc,T_B);
         T=linspace(0.7*Tc,T_inv,q)';

         Pr=zeros(q,1);
         Tr=zeros(q,1);
         
         for i=1:q
             
         [Pr(i), Tr(i), V(i)]=Calculation_JTIC_PR_EoS_For_Pure(alpha_Function,p,Pc,Tc,T(i),w);
%          P(i)=Pc*Pr(i)*1e-05;


         end
         
         
         Tr(q)=T_inv/Tc;Pr(q)=1e-07;
         M=[Tr Pr];

         if q>100
             
         plot(Pr,Tr);xlabel('Pr');ylabel('Tr');
         
         else
             plot(Pr,Tr,'+');xlabel('Pr');ylabel('Tr');
             
         end

         Vc=94;
         vr=V./Vc;Vr=vr';
         Mat=[Tr Pr Vr];
         
% T=Tc.*tr;
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
% [Pr(j,i), Tr(j,i)]=Calculation_JTIC_PR_EoS_For_Pure(alpha_Function,p,Pc,Tc,T(j),w);
% 
% 
%  end
% end
% 
% plot(tr,Pr(:,1),'k',tr,Pr(:,2),'r',tr,Pr(:,3),'b',tr,Pr(:,4),'g')
% legend({'Original','Coquelet','Haghtalab','Saffari'});title(M)
%     
 