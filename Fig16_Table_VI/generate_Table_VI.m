clear variables;
close all;

% Loading data from the file 'Notreatment_max_rho4_10_N=10000' 
load('Notreatment_max_rho4_10_N=10000');

Tnt=Tsurv';

clearvars -except Tnt;

% Loading data from the file 'CAR-T_trial_max_rho4_10_N=10000_v=0_5_10_9'
load('CAR-T_trial_max_rho4_10_N=10000_v=0_5_10_9');


Tsc=Tsurv';

clearvars -except Tnt Tsc;


% Loading data from the file 'TMZ_only_Tsurv_max_r1_0_025_rho4_10_N=10000'
load('TMZ_only_Tsurv_max_r1_0_025_rho4_10_N=10000');

Tst=Tsurv';



clearvars -except Tnt Tsc Tst;


% Loading data from the file 'TMZ_CAR-T_Tsurv_max_r1_0_025_rho4_10_v_5_N=10000'
load('TMZ_CAR-T_Tsurv_max_r1_0_025_rho4_10_v_5_N=10000');

Tstc=Tsurv';

clearvars -except Tnt Tsc Tst Tstc;

% Loading data from the file 'TMZ_CAR-T_TMZ_CAR-T_Tsurv_max_r1_0_025_rho4_10_N=10000' 
load('TMZ_CAR-T_TMZ_CAR-T_Tsurv_max_r1_0_025_rho4_10_N=10000');

Tstctc=Tsurv';

clearvars -except Tnt Tsc Tst Tstc Tstctc;

% Loading data from the file 'CAR-T_TMZ_Tsurv_max_r1_0_025_rho4_10_N=10000'
load('CAR-T_TMZ_Tsurv_max_r1_0_025_rho4_10_N=10000');

Tsct=Tsurv';

clearvars -except Tnt Tsc Tst Tstc Tstctc Tsct;

% Loading data from the file 'TMZ_CAR-T_TMZ_Tsurv_max_r1_0_025_rho4_10_N=10000'
load('TMZ_CAR-T_TMZ_Tsurv_max_r1_0_025_rho4_10_N=10000');

Tstct=Tsurv';

clearvars -except Tnt Tsc Tst Tstc Tstctc Tsct Tstct;

% Loading data from the file 'CAR-T_TMZ_CAR-T_TMZ_surv_max_r1_0_025_rho4_10_N=10000'
load('CAR-T_TMZ_CAR-T_TMZ_surv_max_r1_0_025_rho4_10_N=10000');

Tsctct=Tsurv';

clearvars -except Tnt Tsc Tst Tstc Tstctc Tsct Tstct Tsctct;

% Loading data from the file 'CAR-T_TMZ_CAR-T_surv_max_r1_0_025_rho4_10_v_5_N=10000'
load('CAR-T_TMZ_CAR-T_surv_max_r1_0_025_rho4_10_v_5_N=10000');

Tsctc=Tsurv';

clearvars -except Tnt Tsc Tst Tstc Tstctc Tsct Tstct Tsctct Tsctc;





Times=[Tnt'; Tsc'; Tst'; Tstc'; Tstctc'; Tsct'; Tstct'; Tsctct'; Tsctc'];
% Median values for "Times" variable
for i=1:9
 Ts(i)=round(median(Times(i,:)));
end

clearvars -except Ts;

% Loading data from the file 'Notreatment_max_rho4_10_N=10000'
load('Notreatment_max_rho4_10_N=10000');

Tnt=Tsurv';

clearvars -except Tnt Ts;

% Loading data from the file 'CAR-T_trial_max_rho4_10_v_10_N=10000_v=0_5_10_9'
load('CAR-T_trial_max_rho4_10_v_10_N=10000_v=0_5_10_9');


Tsc=Tsurv';

clearvars -except Tnt Tsc Ts;


% Loading data from the file 'TMZ_only_Tsurv_max_r1_0_025_rho4_10_N=10000' 
load('TMZ_only_Tsurv_max_r1_0_025_rho4_10_N=10000');

Tst=Tsurv';



clearvars -except Tnt Tsc Tst Ts;


% Loading data from the file 'TMZ_CAR-T_Tsurv_max_r1_0_025_rho4_10_v_10_N=10000'
load('TMZ_CAR-T_Tsurv_max_r1_0_025_rho4_10_v_10_N=10000');

Tstc=Tsurv';

clearvars -except Tnt Tsc Tst Tstc Ts;

% Loading data from the file 'TMZ_CAR-T_TMZ_CAR-T_Tsurv_max_r1_0_025_rho4_10_v_10_N=10000'
load('TMZ_CAR-T_TMZ_CAR-T_Tsurv_max_r1_0_025_rho4_10_v_10_N=10000');

Tstctc=Tsurv';

clearvars -except Tnt Tsc Tst Tstc Tstctc Ts;

% Loading data from the file 'CAR-T_TMZ_Tsurv_max_r1_0_025_rho4_10_v_10_N=10000'
load('CAR-T_TMZ_Tsurv_max_r1_0_025_rho4_10_v_10_N=10000');

Tsct=Tsurv';

clearvars -except Tnt Tsc Tst Tstc Tstctc Tsct Ts;

% Loading data from the file 'TMZ_CAR-T_TMZ_Tsurv_max_r1_0_025_rho4_10_v_10_N=10000'
load('TMZ_CAR-T_TMZ_Tsurv_max_r1_0_025_rho4_10_v_10_N=10000');

Tstct=Tsurv';

clearvars -except Tnt Tsc Tst Tstc Tstctc Tsct Tstct Ts;

% Loading data from the file 'CAR-T_TMZ_CAR-T_TMZ_surv_max_r1_0_025_rho4_10_v_10_N=10000'
load('CAR-T_TMZ_CAR-T_TMZ_surv_max_r1_0_025_rho4_10_v_10_N=10000');

Tsctct=Tsurv';

clearvars -except Tnt Tsc Tst Tstc Tstctc Tsct Tstct Tsctct Ts;

% Loading data from the file 'CAR-T_TMZ_CAR-T_surv_max_r1_0_025_rho4_10_v_10_N=10000'
load('CAR-T_TMZ_CAR-T_surv_max_r1_0_025_rho4_10_v_10_N=10000');

Tsctc=Tsurv';

clearvars -except Tnt Tsc Tst Tstc Tstctc Tsct Tstct Tsctct Tsctc Ts;


Times=[Tnt'; Tsc'; Tst'; Tstc'; Tstctc'; Tsct'; Tstct'; Tsctct'; Tsctc'];
% Median values for "Times" variable
for i=1:9
 Ts1(i)=round(median(Times(i,:)));
end


Protocols={'NT','2C','10T','10T2C','5T1C5T1C','2C10T','5T2C5T','1C5T1C5T','1C10T1C'};

Table_VI= table(Protocols',Ts',Ts1');