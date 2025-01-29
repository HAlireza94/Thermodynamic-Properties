function [Tc, Pc, w, Zc]=Critical_Properties(M)
switch M
    %The Temperature and Pressure unit is Kelvin and Bar Respectively    
    
    %Alkanes
    
    case 'Methane'
    Tc=190.56;pc=45.99200;w=0.011;Zc=0.2863;
    case 'Ethane'
    Tc=305.3;pc=49;w=0.099;Zc=0.279;
    case 'Propane'
    Tc=369.9;pc=42.5;w=0.153;Zc=0.276;
    case 'Butane' 
    Tc=425.12;pc=37.96;w=0.2;Zc=0.274;
    case 'Pentane'
    Tc=469.7;pc=33.7;w=0.252;Zc=0.270;
    case 'Hexane'
    Tc=507.6;pc=30.25;w=0.301;Zc=0.266;
    case 'Heptane'
    Tc=540.2;pc=27.4;w=0.350;Zc=0.261;
    case 'Octane'
    Tc=568.7;pc=24.9;w=0.400;Zc=0.256;
    case 'Nonane'
    Tc=594.6;pc=22.9;w=0.444;Zc=0.235;
    case 'Decane'
    Tc=617.7;pc=21.1;w=0.492;Zc=0.231;
    case 'undecane'
    Tc=639;pc=19.5;w=0.530;Zc=0.227;
    case 'dodecane'
    Tc=658.00;pc=18.20;w=0.576;Zc=0.225;
    case 'tridecane'
    Tc=675.00;pc=16.8;w=0.617;Zc=0.220;
    case 'tetradecane'
    Tc=693.00;pc=15.70;w=0.643;Zc=0.226;
    case 'pentadecane'
    Tc=708.00;pc=14.80;w=0.686;Zc=0.215;
    case 'hexadecane'
    Tc=723.00;pc=14.00;w=0.717;Zc=0.213;
    case 'heptadecane'
    Tc=734.00;pc=13.40;w=0.770;Zc=0.215;
    case 'octadecane'
    Tc=747.00;pc=12.95;w=0.814;Zc=0.218;
    case 'nonadecane'
    Tc=773.76;pc=12.53;w=0.853;Zc=0.214;
    case 'eicosane'
    Tc=788.59;pc=11.86;w=0.891;Zc=0.210;
    case 'heneicosane'
    Tc=802.92;pc=11.25;w=0.927;Zc=0.206;
    case 'docosane'
    Tc=816.78;pc=10.68;w=0.963;Zc=0.203;
    case 'tricosane'
    Tc=830.21;pc=9.87;w=0.997;Zc=0.199;
    case 'tetracosane'
    Tc=843.24;pc=9.67;w=1.031;Zc=0.196;
    case 'pentacosane'
    Tc=855.90;pc=9.21;w=1.063;Zc=0.193;
    case 'hexacosane'
    Tc=868.21;pc=8.46;w=1.094;Zc=0.189;
    case 'heptacosane'
    Tc=880.19;pc=8.08;w=1.125;Zc=0.186;
    case 'octacosane'
    Tc=891.87;pc=8.02;w=1.155;Zc=0.183;
    case 'nonacosane'
    Tc=903.26;pc=7.68;w=1.184;Zc=0.180;
    case 'triacontane'
    Tc=914.38;pc=7.35;w=1.212;Zc=0.176;
    case 'henteriacontane'
    Tc=925.25;pc=7.98;w=1.240;Zc=0.173;
    case 'dotricontane'
    Tc=935.87;pc=7.65;w=1.267;Zc=0.170;
    case 'tritriacontane'
    Tc=946.27;pc=7.34;w=1.294;Zc=0.167;
    case 'tetratriacontane'
    Tc=956.45;pc=7.05;w=1.320;Zc=0.165;
    case 'pentatriacontane'
    Tc=966.42;pc=6.78;w=1.345;Zc=0.162;
    case 'hexatriacontane'
    Tc=976.19;pc=6.51;w=1.370;Zc=0.159;
    case 'heptatriacontane'
    Tc=985.78;pc=6.27;w=1.395;Zc=0.156;
    case 'octatriaoctane'
    Tc=995.15;pc=6.03;w=1.419;Zc=0.153;
    case 'nonatriacontane'
    Tc=1004.43;pc=5.81;w=1.442;Zc=0.151;
    case 'tetracontane'
    Tc=1013.50;pc=5.60;w=1.466;Zc=0.148;
    case 'hentetracontane'
    Tc=1022.42;pc=5.41;w=1.488;Zc=0.146;
    case 'dotetracontane'
    Tc=1031.18;pc=5.22;w=1.511;Zc=0.143;
    case 'tritetracontane'
    Tc=1039.80;pc=5.03;w=1.533;Zc=0.141;
    case 'tetratetracontane'
    Tc=1048.29;pc=4.86;w=1.555;Zc=0.138;
    case 'pentatetracontane'
    Tc=1056.64;pc=4.70;w=1.576;Zc=0.136;
    case 'hexatetracontane'
    Tc=1064.86;pc=4.54;w=1.597;Zc=0.134;
    case 'heptatetracontane'
    Tc=1072.96;pc=4.37;w=1.618;Zc=0.131;
    case 'octatetracontane'
    Tc=1080.93;pc=4.23;w=1.638;Zc=0.129;
    case 'nonatetracontane'
    Tc=1088.80;pc=4.09;w=1.658;Zc=0.127;
    case 'pentacontane'
    Tc=1096.55;pc=3.95;w=1.678;Zc=0.125;
    case 'hexacontane'
    Tc=1168.68;pc=3.05;w=1.862;Zc=0.111;
    case 'heptacontane'
    Tc=1232.84;pc=2.53;w=2.026;Zc=0.104;
    case 'octacontane'
    Tc=1290.89;pc=2.13;w=2.174;Zc=0.098;
    case 'nonacontane'
    Tc=1344.06;pc=1.81;w=2.310;Zc=0.092;
    case 'hectane'
    Tc=1393.24;pc=1.56;w=2.436;Zc=0.086;
        
    %Cycloalkanes
    
    
    case 'Cyclopropane'
    Tc=397.91;pc=54.95;w=0.127;Zc=0.270;    
    case 'Cyclobutane'
    Tc=459.93;pc=50.61;w=0.189;Zc=0.278;
    case 'Cyclopentane'
    Tc=511.7;pc=45.1;w=0.196;Zc=0.276;
    case 'Cyclohexane'
    Tc=554;pc=40.7;w=0.212;Zc=0.273;
    case 'Cycloheptane'
    Tc=604.20;pc=38.20;w=0.241;Zc=0.268;
    case 'Cyclooctane'
    Tc=647.20;pc=35.60;w=0.252;Zc=0.271;
    case 'Cyclononane'
    Tc=682.00;pc=33.44;w=0.281;Zc=0.270;
    case 'Cyclodecane'
    Tc=709.00;pc=30.40;w=0.292;Zc=0.261;
    case 'Cycloundecane'
    Tc=728.78;pc=27.37;w=0.476;Zc=0.246;
    case 'Cyclododecane'
    Tc=747.23;pc=25.35;w=0.531;Zc=0.242;
    case 'Cyclotridecane'
    Tc=764.06;pc=23.57;w=0.583;Zc=0.238;
    case 'Cyclotetradecane'
    Tc=777.91;pc=21.94;w=0.632;Zc=0.234;
    case 'Cyclopentadecane'
    Tc=791.72;pc=20.51;w=0.679;Zc=0.230;
    case 'Cyclohexadecane'
    Tc=804.08;pc=19.22;w=0.724;Zc=0.226;
    case 'Cycloheptadecane'
    Tc=815.06;pc=18.04;w=0.767;Zc=0.222;
    case 'Cyclooctadecane'
    Tc=824.69;pc=16.96;w=0.808;Zc=0.218;
    case 'Cyclononadecane'
    Tc=834.39;pc=16.00;w=0.847;Zc=0.214;
    case 'Cycloeicosane'
    Tc=842.81;pc=15.10;w=0.885;Zc=0.211;
    case 'Methyl Cyclopentane'
    Tc=533.2;pc=38;w=0.231;Zc=0.272;
    case 'Ethyl Cyclopentane'
    Tc=569.50;pc=33.9;w=0.271;Zc=0.269;
    case 'Methyl cyclohexane'
    Tc=572.3;pc=34.7;w=0.236;Zc=0.27;
    case 'Ethyl cyclohexane'
    Tc=609.00;pc=30;w=0.243;Zc=0.260;


%Aromatics


    case 'Benzene'
    Tc=562;pc=48.9;w=0.212;Zc=0.268;
    case 'Toluene'
    Tc=593;pc=41;w=0.263;Zc=0.264;
    case 'Ethyl benzene'
    Tc=617.9;pc=37.3;w=0.302;Zc=0.263;
    case 'n-Propyl benzene' %C6H6
    Tc=638.6;pc=32.1;w=0.3487;Zc=0.265;
    case 'n-Butyl benzene'
    Tc=660.5;pc=28.9;w=0.394;Zc=0.262;
    case 'm-Xylene'
    Tc=618;pc=35.4;w=0.325;Zc=0.259;
    
    
    %Alkenes
    
    
    case 'Ethene'
    Tc=282.5;pc=50.42;w=0.087;Zc=0.281;
    case 'Propene'
    Tc=365;pc=46.1;w=0.144;Zc=0.281;  
    case '1-Butene'
    Tc=419.59;pc=40.2;w=0.187;Zc=0.276;
    case '1-Pentene'
    Tc=464.8;pc=35.6;w=0.236;Zc=0.275;
    case '1-Hexene' %C2H5OH
    Tc=504;pc=32.1;w=0.285;Zc=0.272;
    case '1-Heptene'
    Tc=537.30;pc=29.20;w=0.344;Zc=0.267;
    case '1-Octene'
    Tc=567.00;pc=26.80;w=0.392;Zc=0.266;
    case '1-Nonene'
    Tc=593.48;pc=23.30;w=0.374;Zc=0.248;
    case '1-Decene'
    Tc=618.16;pc=22.20;w=0.417;Zc=0.252;
    case '1-Eicosene'
    Tc=801.49;pc=11.07;w=0.885;Zc=0.211;
    case '1-Triacontene'
    Tc=927.80;pc=6.88;w=1.208;Zc=0.177;
    case '1-Tetracontene'
    Tc=1027.21;pc=5.23;w=1.462;Zc=0.149;
    case '1-Pentacontene'
    Tc=1110.43;pc=3.69;w=1.675;Zc=0.125;
    case '1-Hexacontene'
    Tc=1182.70;pc=2.84;w=1.860;Zc=0.111;
    case '1-Heptacontene'
    Tc=1246.96;pc=2.35;w=2.024;Zc=0.104;
    case '1-Octacontene'
    Tc=1305.08;pc=1.98;w=2.172;Zc=0.098;
    
    
    %Esters
    
    
%     case 'Methyl methanoate'
%     Tc=487.20;Pc=59.98;w=0.254;Zc=0.255;
%     case 'Methyl ethanoate'
%     Tc=506.55;Pc=47.50;w=0.331;Zc=0.257;
%     case 'Methyl propanoate' 
%     Tc=530.6;Pc=40.04;w=0.347;Zc=0.256;
%     case 'Methyl butanoate'
%     Tc=554.5;Pc=34.73;w=0.378;Zc=0.256;
%     case 'Methyl pentanoate'
%     Tc=564.61;Pc=29.32;w=0.615;Zc=0.250;
%     case 'Methyl hexanoate'
%     Tc=586.34;Pc=26.08;w=0.682;Zc=0.245;
%     case 'Methyl heptanoate'
%     Tc=606.49;Pc=23.35;w=0.744;Zc=0.241;
%     case 'Methyl octanoate'
%     Tc=625.31;Pc=21.05;w=0.802;Zc=0.237;
%     case 'Methyl nonanoate'
%     Tc=643.00;Pc=19.27;w=0.856;Zc=0.233;
%     case 'Methyl decanoate'
%     Tc=659.70;Pc=17.52;w=0.907;Zc=0.229;
%     case 'Methyl undecanoate'
%     Tc=675.55;Pc=16.04;w=0.956;Zc=0.225;
%     case 'Methyl dodecanoate'
%     Tc=712.15;Pc=15.35;w=1.002;Zc=0.221;
%     case 'Methyl tridecanoate'
%     Tc=705.04;Pc=13.63;w=1.046;Zc=0.217;
%     case 'Methyl tetradecanoate'
%     Tc=718.85;Pc=12.88;w=1.088;Zc=0.213;
%     case 'Methyl pentadecanoate'
%     Tc=732.10;Pc=11.96;w=1.129;Zc=0.210;
%     case 'Methyl hexadecanoate'
%     Tc=744.86;Pc=11.27;w=1.168;Zc=0.206;
%     case 'Methyl heptadecanoate'
%     Tc=757.16;Pc=10.62;w=1.206;Zc=0.202;
%     case 'Methyl octadecanoate'
%     Tc=769.05;Pc=9.77;w=1.242;Zc=0.199;
%     case 'Methyl nonadecanoate'
%     Tc=780.55;Pc=9.53;w=1.278;Zc=0.196;
%     case 'Methyl eicosanoate'
%     Tc=791.69;Pc=9.04;w=1.312;Zc=0.192;
%     case 'Methyl heneicosanoate'
%     Tc=802.50;Pc=8.58;w=1.345;Zc=0.189;
%     case 'Methyl docosanoate'
%     Tc=813.01;Pc=8.15;w=1.377;Zc=0.186;
%     case 'Methyl tetracosanoate'
%     Tc=833.18;Pc=7.37;w=1.439;Zc=0.179;
%     case 'Methyl hexacosanoate'
%     Tc=852.34;Pc=6.68;w=1.498;Zc=0.173;
    
%     case 'Methyl oleate'
%     Tc=764.00;Pc=11.95;w=1.237;Zc=0.199;
    
%     case 'Ethyl methanoate' 
%     Tc=508.40;Pc=47.42;w=0.285;Zc=0.257;
%     case 'Ethyl ethanoate'
%     Tc=523.30;Pc=38.8;w=0.366;Zc=0.255;
%     case 'Ethyl propanoate'
%     Tc=546;Pc=33.62;w=0.394;Zc=0.256;
%     case 'Ethyl butanoate' 
%     Tc=571;Pc=29.21;w=0.615;Zc=0.250;
%     case 'Ethyl pentanoate'
%     Tc=586.34;Pc=25.67;w=0.682;Zc=0.245;
%     case 'Ethyl hexanoate'
%     Tc=606.49;Pc=23.08;w=0.744;Zc=0.241;
%     case 'Ethyl heptanoate'
%     Tc=625.31;Pc=21.17;w=0.802;Zc=0.237;
%     case 'Ethyl octanoate'
%     Tc=643.00;Pc=18.90;w=0.856;Zc=0.233;
%     case 'Ethyl nonanoate'
%     Tc=659.70;Pc=17.34;w=0.907;Zc=0.229;
%     case 'Ethyl decanoate'
%     Tc=675.55;Pc=15.98;w=0.956;Zc=0.225;
%     case 'Ethyl undecanoate'
%     Tc=690.63;Pc=14.77;w=1.002;Zc=0.221;
%     case 'Ethyl dodecanoate'
%     Tc=705.04;Pc=13.70;w=1.046;Zc=0.217;
%     case 'Ethyl tridecanoate'
%     Tc=718.85;Pc=12.63;w=1.088;Zc=0.213;
%     case 'Ethyl tetradecanoate'
%     Tc=732.10;Pc=11.88;w=1.129;Zc=0.210;
%     case 'Ethyl hexadecanoate'
%     Tc=757.16;Pc=10.58;w=1.206;Zc=0.202;
%     case 'Ethyl heptadecanoate'
%     Tc=769.05;Pc=9.79;w=1.242;Zc=0.199;
%     case 'Ethyl octadecanoate'
%     Tc=780.55;Pc=9.53;w=1.278;Zc=0.196;
%     case 'Ethyl eicosanoate'
%     Tc=802.50;Pc=8.58;w=1.345;Zc=0.189;
%     case 'Ethyl docosanoate'
%     Tc=823.23;Pc=7.75;w=1.409;Zc=0.182;
    
%       case 'Propyl methanoate'
%     Tc=538.00;Pc=40.63;w=0.318;Zc=0.259;
%     case 'Butyl methanoate'
%         Tc=540.95;Pc=33.44;w=0.543;Zc=0.254;
%     case 'Pentyl methanoate'
%         Tc=576.00;Pc=31.25;w=0.528;Zc=0.262;
%     case 'Hexyl methanoate'
%         Tc=586.34;Pc=25.74;w=0.682;Zc=0.245;
%     case 'Heptyl methanoate'
%         Tc=606.49;Pc=23.39;w=0.744;Zc=0.241;
%     case 'Octyl methanoate'
%     Tc=625.31;Pc=20.89;w=0.802;Zc=0.237;
%     case 'Nonyl methanoate'
%         Tc=643.00;Pc=19.19;w=0.856;Zc=0.233;
%     case 'Decyl methanoate'
%     Tc=659.70;Pc=17.51;w=0.907;Zc=0.229;
%     case 'Undecyl methanoate'
%         Tc=675.55;Pc=16.04;w=0.956;Zc=0.225;
%     case 'Dodecyl methanoate'
%         Tc=690.63;Pc=14.76;w=1.002;Zc=0.221;
%     case 'Tridecyl methanoate'
%         Tc=705.04;Pc=13.63;w=1.046;Zc=0.217;
%     case 'Tetradecyl methanoate'
%         Tc=718.85;Pc=12.54;w=1.088;Zc=0.213;
%     case 'Pentadecyl methanoate'
%         Tc=732.10;Pc=11.86;w=1.129;Zc=0.210;
%     case 'Hexadecyl methanoate'
%         Tc=744.86;Pc=11.22;w=1.168;Zc=0.206;
%     case 'Heptadecyl methanoate'
%     Tc=757.16;Pc=10.62;w=1.206;Zc=0.202;
%     case 'Octadecyl methanoate'
%     Tc=769.05;Pc=10.06;w=1.242;Zc=0.199;
%     case 'Nonadecyl methanoate'
%         Tc=780.55;Pc=9.53;w=1.278;Zc=0.196;
%     case 'Eicosyl methanoate'
%     Tc=791.69;Pc=9.04;w=1.312;Zc=0.192;
    
    
%      case 'Propyl ethanoate'
%         Tc=549.73;Pc=33.60;w=0.389;Zc=0.254;
%     case 'Butyl ethanoate'
%         Tc=575.40;Pc=30.90;w=0.439;Zc=0.261;
%     case 'Pentyl ethanoate'
%         Tc=599.90;Pc=27.70;w=0.448;Zc=0.259;
%     case 'Hexyl ethanoate'
%         Tc=606.49;Pc=22.92;w=0.744;Zc=0.241;
%     case 'Heptyl ethanoate'
%         Tc=625.31;Pc=20.70;w=0.802;Zc=0.237;
%     case 'Octyl ethanoate'
%         Tc=643.00;Pc=18.83;w=0.856;Zc=0.233;
%     case 'Nonyl ethanoate'
%         Tc=659.70;Pc=17.30;w=0.907;Zc=0.229;
%     case 'Decyl ethanoate'
%         Tc=675.55;Pc=16.04;w=0.956;Zc=0.225;
%     case 'Undecyl ethanoate'
%         Tc=690.63;Pc=14.76;w=1.002;Zc=0.221;
%     case 'Dodecyl ethanoate'
%         Tc=705.04;Pc=13.63;w=1.046;Zc=0.217;
%     case 'Tridecyl ethanoate'
%         Tc=718.85;Pc=12.54;w=1.088;Zc=0.213;
%     case 'Tetradecyl ethanoate'
%         Tc=732.10;Pc=11.86;w=1.129;Zc=0.210;
%     case 'Pentadecyl ethanoate'
%         Tc=744.86;Pc=11.22;w=1.168;Zc=0.206;
%     case 'Hexadecyl ethanoate'
%         Tc=757.16;Pc=10.62;w=1.206;Zc=0.202;
%     case 'Heptadecyl ethanoate'
%         Tc=769.05;Pc=10.06;w=1.242;Zc=0.199;
%     case 'Octadecyl ethanoate'
%         Tc=780.55;Pc=9.53;w=1.278;Zc=0.196;
%     case 'Nonadecyl ethanoate'
%         Tc=791.69;Pc=9.04;w=1.312;Zc=0.192;
%     case 'Eicosyl ethanoate'
%         Tc=802.5;Pc=8.58;w=1.345;Zc=0.189;
        
%     case 'Propyl propanoate'
%         Tc=568.60;Pc=30.60;w=0.449;Zc=0.262;
%     case 'Butyl propanoate'
%         Tc=594.60;Pc=25.99;w=0.682;Zc=0.245;
%     case 'Pentyl propanoate'
%         Tc=606.49;Pc=22.92;w=0.744;Zc=0.241;
%     case 'Hexyl propanoate'
%         Tc=625.31;Pc=20.79;w=0.802;Zc=0.237;
%     case 'Heptyl propanoate'
%         Tc=643.00;Pc=18.84;w=0.856;Zc=0.233;
%     case 'Octyl propanoate'
%         Tc=659.70;Pc=17.45;w=0.907;Zc=0.229;
%     case 'Nonyl propanoate'
%         Tc=675.55;Pc=16.04;w=0.956;Zc=0.225;
%     case 'Decyl propanoate'
%         Tc=690.63;Pc=14.77;w=1.002;Zc=0.221;
%     case 'Undecyl propanoate'
%         Tc=705.04;Pc=13.63;w=1.046;Zc=0.217;
%     case 'Dodecyl propanoate'
%         Tc=718.85;Pc=12.54;w=1.088;Zc=0.213;
%     case 'Tridecyl propanoate'
%         Tc=732.10;Pc=11.86;w=1.129;Zc=0.210;
%     case 'Tetradecyl propanoate'
%         Tc=744.86;Pc=11.22;w=1.168;Zc=0.206;
%     case 'Pentadecyl propanoate'
%         Tc=757.16;Pc=10.62;w=1.206;Zc=0.202;
%     case 'Hexadecyl propanoate'
%         Tc=769.05;Pc=10.06;w=1.242;Zc=0.199;
%     case 'Octadecyl propanoate'
%         Tc=791.69;Pc=9.04;w=1.312;Zc=0.192;
%     case 'Nonadecyl propanoate'
%         Tc=802.50;Pc=8.58;w=1.345;Zc=0.189;
%     case 'Eicosyl propanoate'
%         Tc=813.01;Pc=8.15;w=1.377;Zc=0.186;
        
    case 'Propyl butanoate'
        Tc=593.70;pc=27.20;w=0.433;Zc=0.258;
    case 'Butyl butanoate'
        Tc=606.49;pc=22.95;w=0.744;Zc=0.241;
    case 'Pentyl butanoate'
        Tc=625.31;pc=20.69;w=0.802;Zc=0.237;
    case 'Hexyl butanoate'
        Tc=643.00;pc=18.81;w=0.856;Zc=0.233;
    case 'Heptyl butanoate'
        Tc=659.70;pc=17.10;w=0.907;Zc=0.229;
    case 'Octyl butanoate'
        Tc=675.55;pc=15.76;w=0.956;Zc=0.225;
    case 'Nonyl butanoate'
        Tc=690.63;pc=14.67;w=1.002;Zc=0.221;
    case 'Decyl butanoate'
        Tc=705.04;pc=13.70;w=1.046;Zc=0.217;
    case 'Undecyl butanoate'
        Tc=718.85;pc=12.54;w=1.088;Zc=0.213;
    case 'Dodecyl butanoate'
        Tc=732.10;pc=11.86;w=1.129;Zc=0.210;
    case 'Tridecyl butanoate'
        Tc=744.86;pc=11.22;w=1.168;Zc=0.206;
    case 'Tetradecyl butanoate'
        Tc=757.16;pc=10.62;w=1.206;Zc=0.202;
    case 'Pentadecyl butanoate'
        Tc=769.05;pc=10.06;w=1.242;Zc=0.199;
    case 'Hexadecyl butanoate'
        Tc=780.55;pc=9.53;w=1.278;Zc=0.196;
    case 'Heptadecyl butanoate'
        Tc=791.69;pc=9.04;w=1.312;Zc=0.192;
    case 'Octadecyl butanoate'
        Tc=802.50;pc=8.58;w=1.345;Zc=0.189;
    case 'Nonadecyl butanoate'
        Tc=813.01;pc=8.15;w=1.377;Zc=0.186;
    
%     case 'n-Propyl methanoate'
%     Tc=538;Pc=40.63;w=0.318;Zc=0.259;
%     case 'n-Propyl ethanoate'
%     Tc=549.35;Pc=33.34;w=0.38978;Zc=0.254;
%     case 'n-Butyl ethanoate'
%     Tc=575.40;Pc=30.90;w=0.439;Zc=0.261;
%     case 'n-Propyl propanoate'
%     Tc=568.6;Pc=30.6;w=0.449;Zc=0.262;
    
    
%     case 'n-Propyl butanoate'
%     Tc=209.4;Pc=55;w=0.0;
    
% 
%Ketones


    case 'Acetone' 
    Tc=508;pc=48;w=0.304;Zc=0.233;
    case 'Methyl ethyl ketone'
    Tc=535;pc=42;w=0.320;Zc=0.249;
    case 'Methyl n-propyl ketone'
    Tc=561.1;pc=37;w=0.340;Zc=0.238;
    case 'Diethyl ketone'
    Tc=560.9;pc=37.4;w=0.345;Zc=0.269;
    
    
    %Ethers
    
    
    case 'Dimethyl ether'
    Tc=401;pc=54;w=0.200;Zc=0.274;
    case 'Methyl ethyl ether'
    Tc=437.8;pc=44;w=0.222;Zc=0.267;
    case 'Methyl n-propyl ether'
    Tc=476.25;pc=38.01;w=0.272;Zc=0.265;
    case 'Diethyl ether'
    Tc=467;pc=36;w=0.281;Zc=0.263;
    case 'Phenyl ether'
    Tc=766.8;pc=29.58;w=0.848;Zc=0.233;

%Gases


    case 'Nitrogen'
    Tc=126.2;pc=34;w=0.039;Zc=0.294;
    case 'Argon'
    Tc=150.86;pc=48.98;w=0.000;Zc=0.291;
    case 'Carbon monoxide'
    Tc=134.5;pc=35;w=0.066;Zc=0.295;
    case 'Carbon dioxide'
    Tc=304.2;pc=73.8;w=0.239;Zc=0.274;
    case 'Carbon disulfide'
    Tc=552;pc=79;w=0.109;Zc=0.275;
    case 'Sulfur dioxide'
    Tc=430.3;pc=78.8;w=0.256;Zc=0.269;
    
    
    %Haloalkanes
    
    
    case 'Chloromethane'
    Tc=416.2;pc=66.6;w=0.153;Zc=0.268;
    case 'Dichloromethane'
    Tc=510;pc=63.55;w=0.199;Zc=0.261; 
    case 'Trichloromethane'
    Tc=537.00;pc=53.30;w=0.218;Zc=0.293;
    case 'Tetrachloromethane'
    Tc=556.30;pc=45.40;w=0.193;Zc=0.272;
    case'Chloroethane'
    Tc=460.4;pc=52.7;w=0.191;Zc=0.274;
    case'1-Chloropropane'
    Tc=503.15;pc=45.8;w=0.228;Zc=0.291;   
    
    
    %Alkanols
    
    
    case 'Methanol'
    Tc=513;pc=81;w=0.556;Zc=0.224;
    case 'Ethanol'
    Tc=514;pc=63;w=0.644;Zc=0.24;
    case '1-Propanol'
    Tc=536.78;pc=51.68;w=0.62;Zc=0.254; % with Z=0.7 and T=700 be moshkel mikhord
    case '1-Butanol'
    Tc=563.05;pc=44.24;w=0.591;Zc=0.260;
    case '1-Pentanol'
    Tc=588.1;pc=38.97;w=0.573;Zc=0.260;
    case '1-Hexanol'
    Tc=610.3;pc=34.17;w=0.576;Zc=0.261;
    case '1-Heptanol'
    Tc=632.6;pc=30.58;w=0.567;Zc=0.253;
    case '1-Octanol'
    Tc=652.5;pc=27.77;w=0.583;Zc=0.254;
    case '1-Nonanol'
    Tc=670.7;pc=25.28;w=0.600;Zc=0.259;
    case '1-Decanol'
    Tc=687.3;pc=23.15;w=0.622;Zc=0.263;
    case '2-Propanol'
    Tc=508.3;pc=47.62;w=0.665;Zc=0.248;
    case '2-Butanol'
    Tc=536.00;pc=44.2;w=0.578;Zc=0.252;
    case '2-Pentanol'
    Tc=560.30;pc=36.75;w=0.559;Zc=0.260;
    case '2-Hexanol'
    Tc=585.90;pc=33.10;w=0.556;Zc=0.261;
    case '2-Heptanol'
    Tc=608.32;pc=30.21;w=0.557;Zc=0.264;
    case '2-Octanol'
    Tc=629.60;pc=27.54;w=0.584;Zc=0.273;
    case '2-Nonanol'
    Tc=649.60;pc=25.30;w=0.595;Zc=0.269;
    case '2-Decanol'
    Tc=680.37;pc=20.66;w=0.802;Zc=0.237;
    case '2-Undecanol'
    Tc=698.05;pc=19.03;w=0.856;Zc=0.233;
    case '2-Dodecanol'
    Tc=714.75;pc=17.33;w=0.907;Zc=0.229;
    case '2-Tridecanol'
    Tc=730.59;pc=15.99;w=0.956;Zc=0.225;
    case '2-Tetradecanol'
    Tc=745.68;pc=14.43;w=1.002;Zc=0.221;
    case '2-Pentadecanol'
    Tc=760.09;pc=13.85;w=1.046;Zc=0.217;
    case '2-Hexadecanol'
    Tc=773.89;pc=12.92;w=1.089;Zc=0.213;
    case '2-Heptadecanol'
    Tc=787.14;pc=12.19;w=1.129;Zc=0.210;
    case '2-Octadecanol'
    Tc=799.90;pc=11.53;w=1.168;Zc=0.206;
    case '2-Nonadecanol'
    Tc=812.20;pc=10.91;w=1.206;Zc=0.202;
    case '2-Eicosanol'
    Tc=824.08;pc=10.05;w=1.243;Zc=0.199;
    case '2-Methyl-1-propanol'
    Tc=547.70;pc=43.00;w=0.586;Zc=0.258;
    case '2-Methyl-2-propanol'
    Tc=506.2;pc=39.7;w=0.614;Zc=0.260;
    case 'Phenol'
    Tc=694.3;pc=59.3;w=0.438;Zc=0.243;
        
        
        %Amines
        
        
    case 'Methanamine'
    Tc=430.8;pc=76.2;w=0.283;Zc=0.235;
    case 'Ethanamine'
    Tc=456.45;pc=56.12;w=0.30091;Zc=0.285;
    case '1-Propanamine'
    Tc=499.00;pc=47.40;w=0.283;Zc=0.252;
    case '1-Butanamine'
    Tc=531.90;pc=42;w=0.337;Zc=0.258;
    case '1-Pentanamine'
    Tc=557.66;pc=37.28;w=0.364;Zc=0.265;
    case '1-Hexanamine'
    Tc=589.15;pc=32.88;w=0.399;Zc=0.260;
    case '1-Heptanamine'
    Tc=617.70;pc=29.49;w=0.431;Zc=0.256;
    case '1-Octanamine'
    Tc=643.92;pc=26.61;w=0.460;Zc=0.251;
    case '1-Nonanamine'
    Tc=668.22;pc=24.20;w=0.487;Zc=0.247;
    case '1-Decanamine'
    Tc=690.90;pc=22.18;w=0.538;Zc=0.243;
    case 'Undecylamine'
    Tc=712.21;pc=20.41;w=0.592;Zc=0.238;
    case 'Dodecylamine'
    Tc=732.33;pc=18.95;w=0.644;Zc=0.234;
    case 'Tridecylamine'
    Tc=751.41;pc=17.59;w=0.692;Zc=0.230;
    case 'Tetradecylamine'
    Tc=769.57;pc=16.39;w=0.739;Zc=0.226;
    case 'Pentadecylamine'
    Tc=746.92;pc=14.54;w=0.783;Zc=0.222;
    case 'Hexadecylamine'
    Tc=803.53;pc=14.36;w=0.826;Zc=0.219;
    case 'Heptadecylamine'
    Tc=819.48;pc=14.21;w=0.866;Zc=0.215;
    case 'Octadecylamine'
    Tc=834.82;pc=13.54;w=0.906;Zc=0.211;
    case 'Nonadecylamine'
    Tc=849.62;pc=11.72;w=0.943;Zc=0.207;
    case 'Eicosylamine'
    Tc=863.91;pc=10.98;w=0.980;Zc=0.204;
    case 'Heneicosylamine'
    Tc=877.74;pc=11.03;w=1.015;Zc=0.200;
    case 'Docosylamine'
    Tc=891.14;pc=10.49;w=1.050;Zc=0.197;
    case 'Tricosylamine'
    Tc=904.15;pc=9.99;w=1.083;Zc=0.193;
    case 'Tetracosylamine'
    Tc=916.78;pc=9.52;w=1.115;Zc=0.190;
    case 'Pentacosylamine'
    Tc=929.06;pc=9.07;w=1.146;Zc=0.187;
    case 'Hexacosylamine'
    Tc=941.02;pc=8.65;w=1.177;Zc=0.184;
    case 'Heptacosylamine'
    Tc=952.68;pc=8.25;w=1.207;Zc=0.180;
    case 'Octacosylamine'
    Tc=964.05;pc=7.88;w=1.236;Zc=0.177;
    case 'Nonacosylamine'
    Tc=975.16;pc=7.53;w=1.264;Zc=0.174;
    case 'Triacontylamine'
    Tc=986.01;pc=7.63;w=1.292;Zc=0.171;
    case 'Hentriacontylamine'
    Tc=996.61;pc=7.32;w=1.319;Zc=0.168;
    case 'Dotriacontylamine'
    Tc=1006.99;pc=7.02;w=1.346;Zc=0.165;
    case 'Tetratriacontylamine'
    Tc=1027.11;pc=6.48;w=1.397;Zc=0.160;
    case 'Pentatriacontylamine'
    Tc=1036.88;pc=6.23;w=1.422;Zc=0.157;
    case 'Hexatriacontylamine'
    Tc=1046.45;pc=5.99;w=1.446;Zc=0.154;
    case 'Heptatriacontylamine'
    Tc=1055.85;pc=5.77;w=1.470;Zc=0.152;
    case 'Octatriacontylamine'
    Tc=1065.07;pc=5.55;w=1.494;Zc=0.149;
    case 'Nonatriacontylamine'
    Tc=1074.13;pc=5.35;w=1.517;Zc=0.146;
    case 'Tetracontylamine'
    Tc=1083.04;pc=5.16;w=1.540;Zc=0.144;
    case 'Dimethylamine'
    Tc=437.2;pc=57.12;w=0.227;Zc=0.279;
    case 'Diethylamine'
    Tc=496.6;pc=37.1;w=0.304;Zc=0.270;
    
    
    % Freons
    
    
    case 'R11'
        Tc=471.2;pc=44.08;w=0.189;Zc=0.279;
    case 'R113'
        Tc=487.25;pc=34.1;w=0.252;Zc=0.274;
    case 'R114'
        Tc=410.78;pc=28.17;w=0.268;Zc=0.256;
    case 'R115'
        Tc=353.15;pc=31.57;w=0.251;Zc=0.271;
    case 'R116'
        Tc=292.8;pc=29.8;w=0.249;Zc=0.274;
    case 'R12'
        Tc=384.95;pc=41.25;w=0.18;Zc=0.280;
    case 'R123'
        Tc=456.94;pc=33.76;w=0.282;Zc=0.269;
    case 'R124'
        Tc=395.65;pc=36.6;w=0.288;Zc=0.271;
    case 'R125'
        Tc=339.17;pc=36.2;w=0.305;Zc=0.271;
    case 'R13'
        Tc=301.8;pc=38.7;w=0.172;Zc=0.278;
    case 'R134a'
        Tc=374.18;pc=40.56;w=0.327;Zc=0.259;
    case 'R14'
        Tc=227.5;pc=37.4;w=0.179;Zc=0.277;
    case 'R141b'
        Tc=478.85;pc=43.4;w=0.221;Zc=0.276;
    case 'R142b'
        Tc=410.29;pc=40.41;w=0.231;Zc=0.267;
    case 'R143a'
        Tc=345.88;pc=37.64;w=0.261;Zc=0.255;
    case 'R152a'
        Tc=386.44;pc=45.2;w=0.275;Zc=0.252;
    case 'R21'
        Tc=451.58;pc=51.84;w=0.205;Zc=0.271;
    case 'R218'
        Tc=345.05;pc=26.8;w=0.327;Zc=0.279;
    case 'R22'
        Tc=369.3;pc=49.71;w=0.219;Zc=0.269;
    case 'R227ea'
        Tc=374.83;pc=28.45;w=0.367;Zc=0.250;
    case 'R23'
        Tc=299.01;pc=48.16;w=0.264;Zc=0.256;
    case 'R236ea'
        Tc=412.38;pc=34.2;w=0.369;Zc=0.314;
    case 'R236fa'
        Tc=398.07;pc=32.00;w=0.321;Zc=0.304;
    case 'R245ca'
        Tc=444.57;pc=39.41;w=0.355;Zc=0.324;
    case 'R245fa'
        Tc=427.16;pc=36.51;w=0.33;Zc=0.267;
    case 'R32'
        Tc=351.26;pc=57.84;w=0.227;Zc=0.244;
    case 'R41'
        Tc=317.42;pc=58.75;w=0.198;Zc=0.252;
    case 'RC318'
        Tc=388.37;pc=27.78;w=0.356;Zc=0.279;
        
    case 'R31'
        Tc=411.48;pc=80.04;w=0.155;Zc=0.290;
        
        
%         % Others
%         
%         
    case 'Sulfur hexafluride'
        Tc=318.69;pc=37.6;w=0.2151;Zc=0.282;
    case 'Ammonia'
        Tc=405.4;pc=113;w=0.25;Zc=0.242;
%     case 'Decafluorobutane'
%         Tc=386.35;Pc=23.23;w=0.372;
% %     case 'Dodecafluorobutane'
% %         Tc=694.3;Pc=59.3;w=0.438;
%     case 'Carbonyl sulfide'
%         Tc=378.8;Pc=63.49;w=0.097;
%     case 'Fluorine'
%         Tc=144.12;Pc=51.72;w=0.053;
    case 'Hydrogen sulfide'
        Tc=373.53;pc=89.63;w=0.0942;Zc=0.284;
%     case 'Nitrogen trifluoride'
%         Tc=234.00;Pc=44.61;w=0.12;
%     case 'Propylene'
%         Tc=365.2;Pc=46;w=0.144;
    case 'Propyne'
        Tc=402.4;pc=56.3;w=0.212;Zc=0.276;
%     case 'Isobutane'
%         Tc=407.7;Pc=36.5;w=0.183;
%     case 'Isobutene'
%         Tc=417.9;Pc=39.99;w=0.189;
% %     case 'trans-Butene'
% %         Tc=694.3;Pc=59.3;w=0.438;
%     case 'cis-2-Butene'
%         Tc=435.58;Pc=42.06;w=0.203;
%     case 'Isopentane'
%         Tc=461;Pc=33.8;w=0.227;
%     case 'Neopentane'
%         Tc=433.8;Pc=31.96;w=0.196;
%     case 'Isohexane'
%         Tc=694.3;Pc=59.3;w=0.438;
%     case 'Dodecane'
%         Tc=658.00;Pc=18.2;w=0.576;
% %     case 'Deuterium'
% %         Tc=694.3;Pc=59.3;w=0.438;
% %     case 'Helium'
% %         Tc=694.3;Pc=59.3;w=0.438;
    case 'Hydrogen'
        Tc=33.2;pc=13;w=-0.216;Zc=0.305;
%     case 'Krypton'
%         Tc=209.4;Pc=55.00;w=0.0;
    case 'Neon'
        Tc=44.4;pc=27.6;w=-0.029;Zc=0.3;
%     case 'Nitrous oxide'
%         Tc=309.6;Pc=72.4;w=0.165;
    case 'Oxygen'
        Tc=154.6;pc=50.4;w=0.025;Zc=0.288;
% %     case 'Parahydrogen'
% %         Tc=694.3;Pc=59.3;w=0.438;
    case 'Xenon'
        Tc=289.74;pc=58.4;w=0.0;Zc=0.286;
%         
%         
%         %water
%         
%         
    case 'Water'
        Tc=647.13;pc=220.55;w=0.34449;Zc=0.229;
        
        % Halogen Group
    case 'Fluorine'
    Tc=144.12;pc=51.72;w=0.053;Zc=0.287;
    case 'Chlorine'
    Tc=417.15;pc=77.11;w=0.069;Zc=0.275;
    case 'Bromine'
    Tc=584.15;pc=103.35;w=0.1189;Zc=0.287;
    case 'Iodine'
    Tc=819.15;pc=116.54;w=0.1115;Zc=0.265;
    case 'Astatine'
    Tc=1096.37;pc=1236.30;w=0.6410;Zc=0.3;
    
    % Glycol ether
    
    case '2-Methoxyethanol'
    Tc=564.00;pc=50.84;w=0.387;Zc=0.262;dp=2.129;
    case '2-Ethoxyethanol'
    Tc=588.82;pc=43.12;w=0.475;Zc=0.258;dp=2.081;
    case '2-Propoxyethanol'
    Tc=615.20;pc=36.50;w=0.486;Zc=0.259;dp=2.159;
    case '2-isopropoxyethanol'
    Tc=614.49;pc=36.12;w=0.554;Zc=0.253;
    case '2-Butoxyethanol'
    Tc=633.90;pc=32.70;w=0.521;Zc=0.263;dp=2.081;
    case '2-(2-Methylpropoxy)ethanol'
    Tc=637.84;pc=31.11;w=0.625;Zc=0.249;
    case '2-(Hexyloxy)ethanol'
    Tc=679.27;pc=24.71;w=0.752;Zc=0.240;dp=2.260;
    case '2-(2-Methoxyethoxy)ethanol'
    Tc=630.00;pc=34.45;w=0.635;Zc=0.248;dp=2.749;
    case '2-(2-Ethoxyethoxy)ethanol'
    Tc=656.30;pc=30.53;w=0.7;Zc=0.244;
    case '2-(2-Butoxyethoxy)ethanol'
    Tc=692.30;pc=27.90;w=0.655;Zc=0.260;dp=2.350;
    case '1-Methoxy-2-propanol'
    Tc=588.82;pc=45.20;w=0.475;Zc=0.258;dp=2.401;
    case '1-Propoxy-2-propanol'
    Tc=637.84;pc=31.05;w=0.625;Zc=0.249;dp=1.970;
    case 'propyne'
        Tc=402.40;pc=56.30;w=0.212;Zc=0.276;
        
    case 'furan'
        Tc=490.15;pc=55.0;w=0.202;Zc=0.294;
        
    case 'thiophene'        
          Tc=579.35;pc=56.90;w=0.197;Zc=0.259;
    case 'sulfolane'
        Tc=853.00;pc=58.71;w=0.635;Zc=0.248;
end

Pc=pc*10^5;


end