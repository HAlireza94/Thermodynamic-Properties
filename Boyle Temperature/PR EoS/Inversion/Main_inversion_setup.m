clc
clear
close all

format long
M=input('Enter Name Of Gas :','s');

[Tc, Pc, w]=critical_Properties(M);


disp('==============================================')

EoS=input('name of Equation of State, EoS :','s');
alpha_Function=input('name of Alpha Function, alpha :','s');

disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx')
 
T_B=input('Enter Boyle Temperature, TB=');

disp('<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>')
%  p=linspace(430*10^5,0.001*Pc,100)';
p=1.05*Pc;

 
 
 switch lower(EoS)
     
     case 'pr'
         T_inv=calculation_T_inv(alpha_Function,w,Tc,Pc,T_B);
         T=linspace(0.7*Tc,T_inv,1000)';

         Pr=zeros(1000,1);
         Tr=zeros(1000,1);
         
         for i=1:1000
             
         [Pr(i), Tr(i)]=PR_Equation_of_Satate(alpha_Function,p,Pc,Tc,T(i),w);
         P(i)=Pc*Pr(i)*1e-05;
         end
         
         plot(Pr,Tr);xlabel('Pr');ylabel('Tr');
%          axis([ 2.882e+07 3.797e+06 154.3 618])
         
     case 'srk'
         
     case 'vtpr'
         
     case 'vtsrk'
         
     case 'ptv'
         
 end
 
 
 
 


