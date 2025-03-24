clear variables;
close all;

% Number of points
N=10000;


EventVar=cell(N*2,1);
EventVar(:)={'DECEASED'};

TreatVar=cell(N*2,1);
TreatVar(1:N)={'1C5T1C5T'};
TreatVar(N+1:2*N)={'10T'};

% Loading data from the file 'CAR-T_TMZ_CAR-T_TMZ_surv_max_r1_0_025_r2=r1_rho4_10_v_5_N=10000'
load('CAR-T_TMZ_CAR-T_TMZ_surv_max_r1_0_025_r2=r1_rho4_10_v_5_N=10000');

% Survival time
Tstct=Tsurv';



% Loading data for the file 'TMZ_only_Tsurv_max_r1_0_025_r2=r1_rho4_10_N=10000'
load('TMZ_only_Tsurv_max_r1_0_025_r2=r1_rho4_10_N=10000');

% Survival time
Tstmz=Tsurv';


times=[Tstct;Tstmz];

[p,fh,stats]=MatSurv(times, EventVar,  TreatVar,'TimeMax',2000,'DispP',false,'DispHR',false,'Xlabel','t');








