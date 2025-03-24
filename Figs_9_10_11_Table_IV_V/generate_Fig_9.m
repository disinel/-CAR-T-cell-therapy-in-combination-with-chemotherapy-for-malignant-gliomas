clear variables;
close all;

% Size of the virtual population
N=10000;


EventVar=cell(N*2,1);
EventVar(:)={'DECEASED'};

TreatVar=cell(N*2,1);
TreatVar(1:N)={'2C'};
TreatVar(N+1:2*N)={'NT'};

% Fig 9 upper panel
load('CAR-T_trial_max_rho4_10_v_5_N=10000_v=0_5_10_9');

% Fig 9 bottom panel
%load('CAR-T_trial_max_rho4_10_v_10_N=10000_v=0_5_10_9');

Tsct=Tsurv';

load('Notreatment_max_rho4_10_N=10000');

Tsnt=Tsurv';

times=[Tsct;Tsnt];

[p,fh,stats]=MatSurv(times, EventVar,  TreatVar,'TimeMax',2000,'DispP',false,'DispHR',false,'Xlabel','t');








