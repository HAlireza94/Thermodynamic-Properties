clc
clear
close all
format long

M=input('Enter Name Of Gas :','s');
pin=input('Enter pressure, P=');
p=pin;
[Tc, Pc, w]=Critical_Properties(M);
ssed=experimental_data_for_speed_Sound(M,p);

[A, B, C, D, E]=Heat_Capacity_Parameters(M);


disp('==============================================');

% alpha_function='jub';

disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');

Tx=ssed(:,1);SSx=ssed(:,2);

T=Tx;
P=pin*101325;
a0=input('Enter intial number, a0=');
j=numel(T);s1=numel(P);
Speed_Sound=zeros(numel(T),5);

 for i=1:5
    if i==1
    
        alpha_function='original';
    
    elseif i==2
        
        alpha_function='coquelet';
        
    elseif i==3
        
        alpha_function='haghtalab';
        
    elseif i==4
        
        alpha_function='saffari';
        
    elseif i==5
        alpha_function='jub';

    end
    
    for g=1:j
    
   Speed_Sound(g,i)=Caclulation_Joule_Thomson_Coefficient_VTPR_EoS_For_Pure...
       (alpha_function,Tc,Pc,w,T(g),P,A,B,C,D,E,a0);
   
    end

 end

 plot(T,Speed_Sound(:,1),T,Speed_Sound(:,2),T,Speed_Sound(:,3),T,Speed_Sound(:,4),T,Speed_Sound(:,5),Tx,SSx,'r+');
 
  legend({'Original','Coquelet et al.','Haghtalab et al.','Saffari-Zahedi',...
        'Pina-Martinez et al.'});xlabel('Temperature (K)');ylabel...
        ('Speed of Sound (m/s)');title([M,'Pressure is',p,' atm']);
 
 
 err_o=((abs((Speed_Sound(:,1)-SSx)))./SSx)*100;
 err_c=((abs((Speed_Sound(:,2)-SSx)))./SSx)*100;
 err_h=((abs((Speed_Sound(:,3)-SSx)))./SSx)*100;
 err_s=((abs((Speed_Sound(:,4)-SSx)))./SSx)*100;
 err_j=((abs((Speed_Sound(:,5)-SSx)))./SSx)*100;
 
 total=[mean(err_o) mean(err_c) mean(err_h) mean(err_s) mean(err_j)];disp(total);
 