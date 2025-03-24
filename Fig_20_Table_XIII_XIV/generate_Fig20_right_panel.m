clear variables;
close all;


% Loading data from the file 'TMZ_only_Tsurv_max_r1_0_025_r2=2r1_rho4_10_N=10000' 
load('TMZ_only_Tsurv_max_r1_0_025_r2=2r1_rho4_10_N=10000');

Tst=Tsurv';



clearvars -except Tst;


% Loading data from the file 'TMZ_CAR-T_Tsurv_max_r1_0_025_r2=2r1_rho4_10_v_5_N=10000' 
load('TMZ_CAR-T_Tsurv_max_r1_0_025_r2=2r1_rho4_10_v_5_N=10000');

Tstc=Tsurv';

clearvars -except Tst Tstc;

% Loading data from the file 'TMZ_CAR-T_TMZ_CAR-T_Tsurv_max_r1_0_025_r1=2r2_rho4_10_v_5_N=10000' 
load('TMZ_CAR-T_TMZ_CAR-T_Tsurv_max_r1_0_025_r1=2r2_rho4_10_v_5_N=10000');

Tstctc=Tsurv';

clearvars -except Tnt Tsc Tst Tstc Tstctc;

% Loading data from the file 'CAR-T_TMZ_Tsurv_max_r1_0_025_rho4_r2=2r1_10_v_5_N=10000'
load('CAR-T_TMZ_Tsurv_max_r1_0_025_rho4_r2=2r1_10_v_5_N=10000');

Tsct=Tsurv';

clearvars -except Tnt Tsc Tst Tstc Tstctc Tsct;

% Loading data from the file 'TMZ_CAR-T_TMZ_Tsurv_max_r1_0_025_r2=2r1_rho4_10_v_5_N=10000'
load('TMZ_CAR-T_TMZ_Tsurv_max_r1_0_025_r2=2r1_rho4_10_v_5_N=10000');

Tstct=Tsurv';

clearvars -except Tnt Tsc Tst Tstc Tstctc Tsct Tstct;

% Loading data from the file 'CAR-T_TMZ_CAR-T_TMZ_surv_max_r1_0_025_r2=2r1_rho4_10_v_5_N=10000'
load('CAR-T_TMZ_CAR-T_TMZ_surv_max_r1_0_025_r2=2r1_rho4_10_v_5_N=10000');

Tsctct=Tsurv';

clearvars -except Tst Tstc Tstctc Tstct Tsct Tsctct;

% Loading data from the file 'Notreatment_r2=2r1_max_rho4_10_N=10000'
load('Notreatment_r2=2r1_max_rho4_10_N=10000');

Tnt=Tsurv';

clearvars -except Tst Tstc Tstctc Tstct Tsct Tsctct Tnt;

% Loading data from the file 'CAR-T_trial_r2=2r1_max_rho4_10_v_5_N=10000_'
load('CAR-T_trial_r2=2r1_max_rho4_10_v_5_N=10000_');

Tc=Tsurv';

clearvars -except Tst Tstc Tstctc Tstct Tsct Tsctct Tnt Tc;

% Loading data from the file 'CAR-T_TMZ_CAR-T_surv_r2=2r1_max_r1_0_025_rho4_10_v_5_N=10000'
load('CAR-T_TMZ_CAR-T_surv_r2=2r1_max_r1_0_025_rho4_10_v_5_N=10000');

Tsctc=Tsurv';

clearvars -except Tst Tstc Tstctc Tstct Tsct Tsctct Tnt Tc Tsctc;



Tv=[Tnt; Tst; Tstc; Tstctc;  Tstct; Tsctc; Tsct; Tsctct; Tc];


str(1:10000)="NT"; str(10001:20000)="10T"; str(20001:30000)="10T2C"; str(30001:40000)="5T1C5T1C";  str(40001:50000)="5T2C5T"; 

str(50001:60000)="1C10T1C"; str(60001:70000)="2C10T"; str(70001:80000)="1C5T1C5T";   str(80001:90000)="2C";

str=str';

Str=cellstr(str);

% Fig. 20 Right panel
f=figure();
boxplot(Tv,Str,'Colors','k','Symbol',' ');
ylim([0 1300]);
yticks([0 240 650 1300]);
ylabel('T_{s}');
hold on;
lines1 = findobj(gcf, 'type', 'line'); %,'Tag', 'Box');
set(lines1,'LineWidth',1);
lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
set(lines, 'Color', 'red','LineWidth',1.5);
yline(275,'--black','LineWidth',1.5);
fontsize(f,14,'point');
fontname(f,"Arial");




