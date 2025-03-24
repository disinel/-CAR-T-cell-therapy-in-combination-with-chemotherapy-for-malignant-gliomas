clear variables;
close all;

% Loading data from the file 'Notreatment_max_rho4_10_N=10000' 
load('Notreatment_max_rho4_10_N=10000');

Tnt=Tsurv';

clearvars -except Tnt;

% Loading data from the file 'CAR-T_trial_max_rho4_10_v_10_N=10000_v=0_5_10_9'
load('CAR-T_trial_max_rho4_10_v_10_N=10000_v=0_5_10_9');


Tsc=Tsurv';

clearvars -except Tnt Tsc;


% Loading data from the file 'TMZ_only_Tsurv_max_r1_0_025_rho4_10_N=10000'
load('TMZ_only_Tsurv_max_r1_0_025_rho4_10_N=10000');

Tst=Tsurv';



clearvars -except Tnt Tsc Tst;


% Loading data from the file 'TMZ_CAR-T_Tsurv_max_r1_0_025_rho4_10_v_10_N=10000'
load('TMZ_CAR-T_Tsurv_max_r1_0_025_rho4_10_v_10_N=10000');

Tstc=Tsurv';

clearvars -except Tnt Tsc Tst Tstc;

% Loading data from the file 'TMZ_CAR-T_TMZ_CAR-T_Tsurv_max_r1_0_025_rho4_10_v_10_N=10000' 
load('TMZ_CAR-T_TMZ_CAR-T_Tsurv_max_r1_0_025_rho4_10_v_10_N=10000');

Tstctc=Tsurv';

clearvars -except Tnt Tsc Tst Tstc Tstctc;

% Loading data from the file 'CAR-T_TMZ_Tsurv_max_r1_0_025_rho4_10_v_10_N=10000'
load('CAR-T_TMZ_Tsurv_max_r1_0_025_rho4_10_v_10_N=10000');

Tsct=Tsurv';

clearvars -except Tnt Tsc Tst Tstc Tstctc Tsct;

% Loading data from the file 'CAR-T_TMZ_CAR-T_surv_max_r1_0_025_rho4_10_v_10_N=10000'
load('CAR-T_TMZ_CAR-T_surv_max_r1_0_025_rho4_10_v_10_N=10000');

Tsctc=Tsurv';

clearvars -except Tnt Tsc Tst Tstc Tstctc Tsct Tsctc;

% Loading data from the file 'TMZ_CAR-T_TMZ_Tsurv_max_r1_0_025_rho4_10_v_10_N=10000'
load('TMZ_CAR-T_TMZ_Tsurv_max_r1_0_025_rho4_10_v_10_N=10000');

Tstct=Tsurv';

clearvars -except Tnt Tsc Tst Tstc Tstctc Tsct Tstct Tsctc;

% Loading data from the file 'CAR-T_TMZ_CAR-T_TMZ_surv_max_r1_0_025_rho4_10_v_10_N=10000'
load('CAR-T_TMZ_CAR-T_TMZ_surv_max_r1_0_025_rho4_10_v_10_N=10000');

Tsctct=Tsurv';

clearvars -except Tnt Tsc Tst Tstc Tstctc Tsct Tstct Tsctct Tsctc;




Tv=[Tnt; Tsc; Tst; Tstc; Tstctc; Tsct; Tsctc; Tstct; Tsctct];


str(1:10000)="NT"; str(10001:20000)="2C"; str(20001:30000)="10T";  str(30001:40000)="10T2C"; 

str(40001:50000)="5T1C5T1C";  str(50001:60000)="2C10T"; str(60001:70000)="1C10T1C"; str(70001:80000)="5T2C5T";  str(80001:90000)="1C5T1C5T";  

str=str';

Str=cellstr(str);

% Fig. 16 Bottom panel
figure();
boxplot(Tv,Str,'Colors','k','Symbol',' ');
ylim([0 2700]);
yticks([0 688 1400 2700]);
ylabel('T_{s}');
hold on;
lines1 = findobj(gcf, 'type', 'line');
set(lines1,'LineWidth',1);
yline(688,'--black','LineWidth',1.5);
lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
set(lines, 'Color', 'red','LineWidth',1.5);







