clc
clear
close all

M=input('Enter Name Of Gas :','s');

[Tc, Pc, w]=Critical_Properties(M);

% A=34.942;B=-0.039957;C=0.00019184;D=-1.5303e-07;E=3.9321e-11;Mw=16; % Methane

% A=28.146;B=0.043447;C=0.00018946;D=-1.9082e-07;E=5.3349e-11; % Ethane

% A=31.78;B=0.74489;C=-0.00010945;D=-2.2668e-07;E=9.3458e-11; % Decane

% A=32.083;B=-0.014831;C=0.00024774;D=-2.3766e-07;E=6.8274e-11; % Ethene

% A=31.298;B=0.072449;C=0.00019481;D=-2.1582e-07;E=6.2974e-11; % Propene

% A=29.414;B=-0.0045993;C=0.000013004;D=-5.4759e-09;E=2.9239e-13; % Nitrogen

% A=29.526;B=-8.8999e-03;C=3.8083e-05;D=-3.2629e-08;E=8.8607e-12; % Oxygen
 
% A=27.437;B=0.042315;C=-0.000019555;D=3.9968e-09;E=-2.9872e-13;Mw=44; % CO2

% A=29.556;B=-0.0065807;C=0.00002013;D=-1.2227e-08;E=2.2617e-12; % CO

%  A=20.786;B=0;C=0;D=0;E=0; % Argon

%  A=20.786;B=0;C=0;D=0;E=0; % Xenon

% A=15.278;B=0.19916;C=-0.00016369;D=5.1686e-08;E=-3.182e-12; % Tetrafluoromethane ( R - 14 )

% A=45.579;B=0.49467;C=-0.00040808;D=1.3789e-07;E=-1.1769e-11; % Octafluorocyclobutane (RC-318)

% A=-31.368;B=4.75e-01;C=-3.11e-04;D=8.52e-08;E=-5.05e-12; % Benzene

% A=40.046;B=-0.038287;C=0.00024529;D=-2.1679E-07;E=5.9909E-11; % Methanol

% A=33.573;B=-0.012581;C=0.000088906;D=-7.1783e-08;E=1.8569e-11; % Ammonia

% A=33.933;B=-0.0084186;C=0.000029906;D=-1.7825e-08;E=3.6934e-12;Mw=18; % Water

% A=71.498;B=7.2559*10^-1;C=1.1553*10^-4;D=-4.1196*10^-7;E=1.414*10^-10;Mw=170; %Dodecane

% A=131.75;B=6.7397*10^-1;C=8.7770*10^-4;D=-1.2430*10^-6;E=3.9785*10^-10;Mw=226; %Hexadecane

A=137.73;B=1.0992;C=3.6839*10^-4;D=-8.2058*10^-7;E=2.7259*10^-10;Mw=282; %Hexadecane

disp('==============================================');

% alpha_function=input('name of Alpha Function, alpha :','s');

disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
% % 



T=linspace(Tc+1,1000,300)'; % K
% Mx=[];
% Tx=[];

% cvx=Mx(:,1);
% cpx=Mx(:,2);ssx=Mx(:,3);jtx=Mx(:,4);
% jtx=Mx(:,1);Tx=Mx(:,2);
a0=input('Enter intial number, a0=');
j=numel(T);

P=360*101325;
% P=100*10^5;
s1=numel(P);
% J_T=zeros(j,s1);
% 
% for s=1:s1
%     
% for i=1:j
%     
%    J_T(i,s)=Caclulation_Joule_Thomson_Coefficient_PR_EoS_For_Pure...
%        (alpha_function,Tc,Pc,w,T(i),P(s),A,B,C,D,E,a0);
% 
% end
% 
% end
% 
% jt=J_T.*101325;
% 
% if numel(T)==300 || numel(T)==1000
% % plot(T,jt)
%     plot(T,jt(:,1),T,jt(:,2),T,jt(:,3))
%     legend({'300','360','420'})
% else
% % plot(T,jt)
%     plot(T,jt(:,1),T,jt(:,2),T,jt(:,3))
%     legend({'300','360','420'})
%     
% end

% J_T=zeros(j,4);
% 
% for i=1:4
%     if i==1
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
%     for g=1:j
%     
%    J_T(g,i)=Caclulation_Joule_Thomson_Coefficient_PR_EoS_For_Pure...
%        (alpha_function,Tc,Pc,w,T(g),P,A,B,C,D,E,a0);
%    
%     end
% 
%     
% end
% 
% jt=J_T.*101325;
% 
% if numel(T)==300 || numel(T)==1000
%     plot(T,jt(:,1),T,jt(:,2),T,jt(:,3),T,jt(:,4))
%     legend({'Original','Coquelet','Haghtalab','Saffari'})
% else 
%     plot(T,jt(:,1),'k+',T,jt(:,2),'g*',T,jt(:,3),'o',T,jt(:,4),'r>')
%     legend({'Original','Coquelet','Haghtalab','Saffari'})
%     
% end


J_T=zeros(j,5,s1);
Cp=zeros(j,5,s1);
Cv=zeros(j,5,s1);
Speed_Sound=zeros(j,5,s1);

for q=1:s1
    
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
    
   [J_T(g,i,q), Speed_Sound(g,i,q), Cp(g,i,q), Cv(g,i,q)]=Caclulation_Joule_Thomson_Coefficient_PR_EoS_For_Pure...
       (alpha_function,Tc,Pc,w,T(g),P(q),A,B,C,D,E,a0,Mw);
   
    end

    
end

end
jt=J_T.*101325;

if numel(T)==300 || numel(T)==1000
  figure (1)
    plot(T,jt(:,1,1),T,jt(:,2,1),T,jt(:,3,1),T,jt(:,4,1),T,jt(:,5,1))
    legend({'Original','Coquelet et al.','Haghtalab et al.','Saffari-Zahedi',...
        'Pina-Martinez et al.'});xlabel('Temperature (K)');ylabel...
        ('Jould - Thomson coefficient (K/bar)');title('using PR EoS for Eicosane at P= 360 atm')
%     figure (2)
%     plot(T,Cp(:,1,1),T,Cp(:,2,1),T,Cp(:,3,1),T,Cp(:,4,1),T,Cp(:,5,1),Tx,cpx,'r+')
%     legend({'Original','Coquelet et al.','Haghtalab et al.','Saffari-Zahedi',...
%         'Pina-Martinez et al.','Experimental data'});xlabel('Temperature (K)');ylabel...
%         ('isobaric heat capacity (Cp)');title('using VTPR EoS for Methane at P= 100 bar')
%     figure (3)
%     plot(T,Cv(:,1,1),T,Cv(:,2,1),T,Cv(:,3,1),T,Cv(:,4,1),T,Cv(:,5,1),Tx,cvx,'r+')
%     legend({'Original','Coquelet et al.','Haghtalab et al.','Saffari-Zahedi',...
%         'Pina-Martinez et al.','Experimental data'});xlabel('Temperature (K)');ylabel...
%         ('isochoric heat capacity (Cv)');title('using VTPR EoS for Methane at P= 100 bar') %### I guess it is wrong !!!
%     figure (4)
%     plot(T,Speed_Sound(:,1,1),T,Speed_Sound(:,2,1),T,Speed_Sound(:,3,1),T,Speed_Sound(:,4,1)...
%         ,T,Speed_Sound(:,5,1),Tx,ssx,'r+')
%     legend({'Original','Coquelet et al.','Haghtalab et al.','Saffari-Zahedi',...
%         'Pina-Martinez et al.','Experimental data'});xlabel('Temperature (K)');ylabel...
%         ('Speed of Sound (m/s)');title('using VTPR EoS for Carbon dioxide at P= 100 bar')
%     
    %     figure (2)
%     plot(T,jt(:,1,2),T,jt(:,2,2),T,jt(:,3,2),T,jt(:,4,2))
%     legend({'Original','Coquelet','Haghtalab','Saffari'});title('P=480 atm')
    
%     figure (3)
%     plot(T,jt(:,1,3),T,jt(:,2,3),T,jt(:,3,3),T,jt(:,4,3))
%     legend({'Original','Coquelet','Haghtalab','Saffari'});title('P=420 atm')
    
else 
    
    figure (1)
    plot(T,jt(:,1,1),'k+',T,jt(:,2,1),'g*',T,jt(:,3,1),'bo',T,jt(:,4,1),'r>')
    legend({'Original','Coquelet','Haghtalab','Saffari'});title('P=420 atm')
    
%     figure (2)
%     plot(T,jt(:,1,2),'k+',T,jt(:,2,2),'g*',T,jt(:,3,2),'bo',T,jt(:,4,2),'r>')
%     legend({'Original','Coquelet','Haghtalab','Saffari'});title('P=480 atm')
    
%     figure (3)
%     plot(T,jt(:,1,3),T,jt(:,2,3),T,jt(:,3,3),T,jt(:,4,3))
%     legend({'Original','Coquelet','Haghtalab','Saffari'});title('P=420 atm')
    
end

% M480=jt(:,:,2);
% M420=jt(:,:,1);