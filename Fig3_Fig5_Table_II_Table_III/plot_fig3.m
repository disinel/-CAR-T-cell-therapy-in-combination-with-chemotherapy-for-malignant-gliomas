clear variables;
close all;

% Load the data from computations

load('TMZ_only_Tsurv_max_r1_0_025_rho4_10_N=10000');

Tstmz=Tsurv';

clearvars -except Tstmz;

load('Notreatment_max_rho4_10_N=10000');

Tsnt=Tsurv';

clearvars -except Tstmz Tsnt N;

n_patinets=N;

% Preparing data for MatSurv format


EventVar=cell(n_patinets*2,1);
EventVar(:)={'DECEASED'};

TreatVar=cell(n_patinets*2,1);
TreatVar(1:n_patinets)={'10T'};
TreatVar(n_patinets+1:2*n_patinets)={'NT'};

times=[Tstmz;Tsnt];

% Plotting K-M curves with MatSurv
[p,fh,stats]=MatSurv(times, EventVar,  TreatVar,'TimeMax',2000,'DispP',false,'DispHR',false,'Xlabel','t');








