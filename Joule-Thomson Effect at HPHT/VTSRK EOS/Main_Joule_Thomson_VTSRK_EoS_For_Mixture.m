clc
clear
close all

yi=[5.651
0.284
83.3482
7.526
2.009
0.305
0.52
0.12
0.144
0.068
0.0138
0.011];Yi=yi./sum(yi);

Kij=[0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0];

fg=[29.414	-0.0045993	0.000013004	-5.4759E-09	2.9239E-13
27.437	0.042315	-0.000019555	3.9968E-09	-2.9872E-13
34.942	-0.039957	0.00019184	-1.5303E-07	3.9321E-11
28.146	0.043447	0.00018946	-1.9082E-07	5.3349E-11
28.277	0.116	0.00019597	-2.3271E-07	6.8669E-11
6.772	3.41E-01	-1.03E-04	-3.68E-08	2.04E-11
20.056	0.28153	-0.000013143	-9.4571E-08	3.4149E-11
-0.881	4.75E-01	-2.48E-04	6.75E-08	-8.53E-12
26.671	0.32324	0.00004282	-1.6639E-07	5.6036E-11
25.924	0.41927	-0.000012491	-1.5916E-07	5.8784E-11
26.984	0.50387	-0.000044748	-1.6835E-07	6.5183E-11
29.053	0.58016	-0.000057103	-1.9548E-07	7.6614E-11];

A=fg(:,1);B=fg(:,2);C=fg(:,3);D=fg(:,4);E=fg(:,5);
Cri_P=[126.2	34.6	0.0377	0.294
304.19	73.82	0.228	0.274
190.56	45.99	0.012	0.286
305.32	48.72	0.1	0.279
369.83	42.48	0.152	0.276
407.8	36.4	0.184	0.278
425.12	37.96	0.2	0.274
460.4	33.8	0.228	0.27
469.7	33.7	0.252	0.27
507.6	30.25	0.301	0.266
540.2	27.4	0.35	0.261
568.7	24.9	0.4	0.256];

Tc=Cri_P(:,1);pc=Cri_P(:,2);Pc=pc.*10^5;
w=Cri_P(:,3);Zc=Cri_P(:,4);

T=linspace(300,600,1000)';

% T=500;


N=numel(Yi);
pcm=zeros(1,N);
tcm=zeros(1,N);
for i=1:N
    
    pcm(i)=Pc(i)*Yi(i);
    tcm(i)=Tc(i)*Yi(i);

end
PcM=sum(pcm);TcM=sum(tcm);
CP=[TcM PcM];


disp('==============================================')

alpha_function=input('name of Alpha Function, alpha :','s');




% T=310.93;
P=270*10^5; % We have to input pressure in bar and then we convert it to pascal. 
% by the way we should pay attention to heat capacity factors at constant
% pressure because we don't have much data and there are limites.
% T=linspace(150,2.5*Tc_m,1000)';

% p=linspace(250,2000,1000)';
j=numel(T);

% P=p.*10^5;
s1=numel(P);
J_T_m=zeros(j,s1);V=zeros(j,s1);

for s=1:s1
    
for i=1:j
    
   [J_T_m(i,s), V(i,s)]=Caclulation_Joule_Thomson_Coefficient_VTSRK_EoS_For_Mixture...
       (alpha_function,T(i),Tc,Pc,w,N,Yi,Kij,A,B,C,D,E,P(s),Zc);

end

end

jt_m=J_T_m.*10^6;
P1=P*10^-6;
M=[T jt_m];
