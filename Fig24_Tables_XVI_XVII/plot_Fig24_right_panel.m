clear variables;
close all;

% Number of points
N=10000;


EventVar=cell(N*2,1);
EventVar(:)={'DECEASED'};

TreatVar=cell(N*2,1);
TreatVar(1:N)={'2C'};
TreatVar(N+1:2*N)={'NT'};

% Loading data from the file 'CAR-T_trial_r2=2r1_max_rho4_10_v_5_N=10000_'
load('CAR-T_trial_r2=2r1_max_rho4_10_v_5_N=10000_');

Tstct=Tsurv';



% Loading data from the file 'Notreatment_r2=2r1_max_rho4_10_N=10000'
load('Notreatment_r2=2r1_max_rho4_10_N=10000');

Tstmz=Tsurv';


times=[Tstct;Tstmz];

[p,fh,stats]=MatSurv(times, EventVar,  TreatVar,'TimeMax',2000,'DispP',false,'DispHR',false,'Xlabel','t');








