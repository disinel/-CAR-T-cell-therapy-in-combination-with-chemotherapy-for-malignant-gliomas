clear variables;
close all;

% Load data from computations

load('CAR-T_trial_max_rho4_5_v_5_N=10000_');

Ts(1)=median(Tsurv);

r4(1)=max(rho4val);

clearvars -except Ts r4;

load('CAR-T_trial_max_rho4_10_N=10000_v=0_5_10_9');

Ts(2)=median(Tsurv);

r4(2)=max(rho4val);

clearvars -except Ts r4;

load('CAR-T_trial_max_rho4_20_v_5_N=10000_v=0_5_10_9');

Ts(3)=median(Tsurv);

r4(3)=max(rho4val);

clearvars -except Ts r4;

load('CAR-T_trial_max_rho4_30_v_5_N=10000_');

Ts(4)=median(Tsurv);

r4(4)=max(rho4val);

clearvars -except Ts r4;

load('CAR-T_trial_max_rho4_40_v_5_N=10000_');

Ts(5)=median(Tsurv);

r4(5)=max(rho4val); 

clearvars -except Ts r4;

load('CAR-T_trial_max_rho4_50_v_5_N=10000_');

Ts(6)=median(Tsurv);

r4(6)=max(rho4val); 

clearvars -except Ts r4;

load('CAR-T_trial_max_rho4_60_v_5_N=10000_');

Ts(7)=median(Tsurv);

r4(7)=max(rho4val); 

clearvars -except Ts r4;

load('CAR-T_trial_max_rho4_70_v_5_N=10000_');

Ts(8)=median(Tsurv);

r4(8)=max(rho4val); 

clearvars -except Ts r4;

load('CAR-T_trial_max_rho4_80_v_5_N=10000_');

Ts(9)=median(Tsurv);

r4(9)=max(rho4val); 

clearvars -except Ts r4;

load('CAR-T_trial_max_rho4_90_v_5_N=10000_');

Ts(10)=median(Tsurv);

r4(10)=max(rho4val); 

clearvars -except Ts r4;

load('CAR-T_trial_max_rho4_100_v_5_N=10000_');

Ts(11)=median(Tsurv);

r4(11)=max(rho4val); 

clearvars -except Ts r4;


load('CAR-T_trial_max_rho4_5_v_10_N=10000_');

Ts1(1)=median(Tsurv);

clearvars -except Ts r4 Ts1;

load('CAR-T_trial_max_rho4_10_v_10_N=10000_');

Ts1(2)=median(Tsurv);

clearvars -except Ts r4 Ts1;

load('CAR-T_trial_max_rho4_20_v_10_N=10000_');

Ts1(3)=median(Tsurv);

clearvars -except Ts r4 Ts1;

load('CAR-T_trial_max_rho4_30_v_10_N=10000_');

Ts1(4)=median(Tsurv);

clearvars -except Ts r4 Ts1;

load('CAR-T_trial_max_rho4_40_v_10_N=10000_');

Ts1(5)=median(Tsurv);

clearvars -except Ts r4 Ts1;

load('CAR-T_trial_max_rho4_50_v_10_N=10000_');

Ts1(6)=median(Tsurv);

clearvars -except Ts r4 Ts1;

load('CAR-T_trial_max_rho4_60_v_10_N=10000_');

Ts1(7)=median(Tsurv);

clearvars -except Ts r4 Ts1;

load('CAR-T_trial_max_rho4_70_v_10_N=10000_');

Ts1(8)=median(Tsurv);

clearvars -except Ts r4 Ts1;

load('CAR-T_trial_max_rho4_80_v_10_N=10000_');

Ts1(9)=median(Tsurv);

clearvars -except Ts r4 Ts1;

load('CAR-T_trial_max_rho4_90_v_10_N=10000_');

Ts1(10)=median(Tsurv);

clearvars -except Ts r4 Ts1;

load('CAR-T_trial_max_rho4_100_v_10_N=10000_');

Ts1(11)=median(Tsurv);

clearvars -except Ts r4 Ts1;


%generating Fig. 6

f=figure();
plot(r4,Ts,'-oblack','Linewidth',1.5);
hold on;
plot(r4,Ts1,'--oblack','Linewidth',1.5);
legend('$2v=10^{9}$','$2v=2\cdot 10^{9}$','Interpreter','latex');
xlim([0.05 1]);
ylim([270 351]);
xticks([0.05 0.5 1]);
yticks([270 310 350]);
ylabel('$\widetilde{T}_{s}^{2C}$','Interpreter','latex');
xlabel('$\max \rho_{4}$','Interpreter','latex');
fontsize(f,14,'point');
fontname(f,"Arial");



