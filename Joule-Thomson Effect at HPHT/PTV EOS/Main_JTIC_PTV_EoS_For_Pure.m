clc
clear
close all

format long

N=input('Enter Name Of Gas :','s');

[Tc, Pc, w, Zc]=Critical_Properties(N);


disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx')
 
T_B=input('Enter Boyle Temperature, TB=');

disp('<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>')
p=1.0*Pc;

 n=input('Enter the number of array, q=');

 
         T_inv=Calculation_T_Inv_PTV_EoS_For_Pure(w,Tc,Pc,T_B,Zc);
         T=linspace(0.7*Tc,T_inv,n)';

% tr=[];
% n=numel(tr);T=Tc.*tr;
         Pr=zeros(n,1);
         Tr=zeros(n,1);
         
         for i=1:n
             
         [Pr(i), Tr(i)]=Calculation_JTIC_PTV_EoS_For_Pure(p,Pc,Tc,T(i),w,Zc);
%          P(i)=Pc*Pr(i)*1e-05;


         end
                  Tr(n)=T_inv/Tc;Pr(n)=1e-07;

         M=[Tr Pr];
         plot(Tr,Pr);xlabel('Pr');ylabel('Tr');

 
 
 


