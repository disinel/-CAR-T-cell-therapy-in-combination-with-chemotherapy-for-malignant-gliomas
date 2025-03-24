clear variables;
close all;

% Loading data from the file 'Notreatment_r2=2r1_max_rho4_10_N=10000'
load('Notreatment_r2=2r1_max_rho4_10_N=10000');

Tnt=Tsurv';

clearvars -except Tnt;

% Loading data from the file 'CAR-T_trial_r2=2r1_max_rho4_10_v_5_N=10000_'
load('CAR-T_trial_r2=2r1_max_rho4_10_v_5_N=10000_');


Tsc=Tsurv';

clearvars -except Tnt Tsc;


% Loading data from the file 'TMZ_only_Tsurv_max_r1_0_025_r2=2r1_rho4_10_N=10000'
load('TMZ_only_Tsurv_max_r1_0_025_r2=2r1_rho4_10_N=10000');

Tst=Tsurv';



clearvars -except Tnt Tsc Tst;


% Loading data from the file 'TMZ_CAR-T_Tsurv_max_r1_0_025_r2=2r1_rho4_10_v_5_N=10000'
load('TMZ_CAR-T_Tsurv_max_r1_0_025_r2=2r1_rho4_10_v_5_N=10000');

Tstc=Tsurv';

clearvars -except Tnt Tsc Tst Tstc;

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

clearvars -except Tnt Tsc Tst Tstc Tstctc Tsct Tstct Tsctct;

% Loading data from the file 'CAR-T_TMZ_CAR-T_surv_r2=2r1_max_r1_0_025_rho4_10_v_5_N=10000'
load('CAR-T_TMZ_CAR-T_surv_r2=2r1_max_r1_0_025_rho4_10_v_5_N=10000');

Tsctc=Tsurv';

clearvars -except Tnt Tsc Tst Tstc Tstctc Tsct Tstct Tsctct Tsctc;

% Survival time
Ts=[round(median(Tnt)), round(median(Tsc)), round(median(Tst)), round(median(Tstc)), round(median(Tstctc)),...
    round(median(Tsct)), round(median(Tstct)), round(median(Tsctct)), round(median(Tsctc))];

clearvars -except Ts;

% Loading data from the file 'Notreatment_r2=2r1_max_rho4_10_N=10000'
load('Notreatment_r2=2r1_max_rho4_10_N=10000');

Tnt=Tsurv';

clearvars -except Ts Tnt;

% Loading data from the file 'CAR-T_trial_r2=2r1_max_rho4_10_v_10_N=10000_'
load('CAR-T_trial_r2=2r1_max_rho4_10_v_10_N=10000_');


Tsc=Tsurv';

clearvars -except Ts Tnt Tsc;


% Loading data from the file 'TMZ_only_Tsurv_max_r1_0_025_r2=2r1_rho4_10_N=10000'
load('TMZ_only_Tsurv_max_r1_0_025_r2=2r1_rho4_10_N=10000');

Tst=Tsurv';



clearvars -except Ts Tnt Tsc Tst;


% Loading data from the file 'TMZ_CAR-T_Tsurv_max_r1_0_025_r2=2r1_rho4_10_v_10_N=10000' 
load('TMZ_CAR-T_Tsurv_max_r1_0_025_r2=2r1_rho4_10_v_10_N=10000');

Tstc=Tsurv';

clearvars -except Ts Tnt Tsc Tst Tstc;

% Loading data from the file 'TMZ_CAR-T_Tsurv_max_r1_0_025_r2=2r1_rho4_10_v_10_N=10000'
load('TMZ_CAR-T_TMZ_CAR-T_Tsurv_max_r1_0_025_r1=2r2_rho4_10_v_10_N=10000');

Tstctc=Tsurv';

clearvars -except Ts Tnt Tsc Tst Tstc Tstctc;

% Loading data from the file 'CAR-T_TMZ_Tsurv_max_r1_0_025_rho4_r2=2r1_10_v_10_N=10000'
load('CAR-T_TMZ_Tsurv_max_r1_0_025_rho4_r2=2r1_10_v_10_N=10000');

Tsct=Tsurv';

clearvars -except Ts Tnt Tsc Tst Tstc Tstctc Tsct;

% Loading data from the file 'TMZ_CAR-T_TMZ_Tsurv_max_r1_0_025_r2=2r1_rho4_10_v_10_N=10000'
load('TMZ_CAR-T_TMZ_Tsurv_max_r1_0_025_r2=2r1_rho4_10_v_10_N=10000');

Tstct=Tsurv';

clearvars -except Ts Tnt Tsc Tst Tstc Tstctc Tsct Tstct;

% Loading data from the file 'CAR-T_TMZ_CAR-T_TMZ_surv_max_r1_0_025_r2=2r1_rho4_10_v_10_N=10000'
load('CAR-T_TMZ_CAR-T_TMZ_surv_max_r1_0_025_r2=2r1_rho4_10_v_10_N=10000');

Tsctct=Tsurv';

clearvars -except Ts Tnt Tsc Tst Tstc Tstctc Tsct Tstct Tsctct;

% Loading data from the file 'CAR-T_TMZ_CAR-T_surv_r2=2r1_max_r1_0_025_rho4_10_v_10_N=10000'
load('CAR-T_TMZ_CAR-T_surv_r2=2r1_max_r1_0_025_rho4_10_v_10_N=10000');

Tsctc=Tsurv';

clearvars -except Ts Tnt Tsc Tst Tstc Tstctc Tsct Tstct Tsctct Tsctc;

% Survival time
Ts1=[round(median(Tnt)), round(median(Tsc)), round(median(Tst)), round(median(Tstc)), round(median(Tstctc)),...
    round(median(Tsct)), round(median(Tstct)), round(median(Tsctct)), round(median(Tsctc))];


% Protocols used in the article
Protocols={'NT','2C','10T','10T2C','5T1C5T1C','2C10T','5T2C5T','1C5T1C5T','1C10T1C'};

Table_X= table(Protocols',Ts',Ts1');