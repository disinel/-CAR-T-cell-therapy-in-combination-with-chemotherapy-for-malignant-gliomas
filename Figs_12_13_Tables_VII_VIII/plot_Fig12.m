clear variables;
close all;

% Size of the virtual population
N=10000;


EventVar=cell(N*2,1);
EventVar(:)={'DECEASED'};

TreatVar=cell(N*2,1);
TreatVar(1:N)={'5T2C5T'};
TreatVar(N+1:2*N)={'10T'};



% Fig. 12 upper panel 
load('TMZ_CAR-T_TMZ_Tsurv_max_r1_0_025_rho4_10_v_5_N=10000');

Tstct=Tsurv';


% Fig. 12 bottom panel
% load('TMZ_CAR-T_TMZ_Tsurv_max_r1_0_025_rho4_10_v_10_N=10000');
% 
% Tstct=Tsurv';


load('TMZ_only_Tsurv_max_r1_0_025_rho4_10_N=10000');

Tstmz=Tsurv';


times=[Tstct;Tstmz];

[p,fh,stats]=MatSurv(times, EventVar,  TreatVar,'TimeMax',2000,'DispP',false,'DispHR',false,'Xlabel','t');








