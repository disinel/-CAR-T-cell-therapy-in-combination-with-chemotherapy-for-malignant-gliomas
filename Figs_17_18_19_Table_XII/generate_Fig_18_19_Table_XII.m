clear variables;
close all;

% Loading data from the file 'TMZ_CAR-T_TMZ_Tsurv_max_r1_0_025_rho4_10_v_5_N=10000' 
load('TMZ_CAR-T_TMZ_Tsurv_max_r1_0_025_rho4_10_v_5_N=10000');

Ttct=Tsurv;

clearvars -except Ttct;

% Loading data from the file 'CAR-T_TMZ_CAR-T_TMZ_surv_max_r1_0_025_rho4_10_v_5_N=10000'
load('CAR-T_TMZ_CAR-T_TMZ_surv_max_r1_0_025_rho4_10_v_5_N=10000');

Tctct=Tsurv;

clearvars -except Ttct Tctct;

% Loading parameters 
load('parameters_N=10000');


Tfrac=Tctct./Ttct;



% Computation of the correlation coefficient presented in Table XII.  
[R,P]=corrcoef(r1val,Tfrac); cf(1)=R(1,2); pv(1)=P(1,2);   [R,P]=corrcoef(alpha1val,Tfrac); cf(2)=R(1,2); pv(2)=P(1,2);  

[R,P]=corrcoef(alpha3val,Tfrac); cf(3)=R(1,2);  pv(3)=P(1,2); [R,P]=corrcoef(epsilon1val,Tfrac); cf(4)=R(1,2);  pv(4)=P(1,2);

[R,P]=corrcoef(rho1val,Tfrac); cf(5)=R(1,2); pv(5)=P(1,2);  [R,P]=corrcoef(rho2val,Tfrac); cf(6)=R(1,2);  pv(6)=P(1,2);

[R,P]=corrcoef(rho3val,Tfrac);  cf(7)=R(1,2); pv(7)=P(1,2); [R,P]=corrcoef(rho4val,Tfrac); cf(8)=R(1,2); pv(8)=P(1,2);

[R,P]=corrcoef(T0val,Tfrac);  cf(9)=R(1,2); pv(9)=P(1,2); [R,P]=corrcoef(delta1val,Tfrac); cf(10)=R(1,2); pv(10)=P(1,2); 

[R,P]=corrcoef(delta2val,Tfrac); cf(11)=R(1,2);  pv(11)=P(1,2);


Vars={'r1','alpha1','alpha3','epsilon1','rho1','rho2','rho3','rho4','T0','delta1','delta2'};

Table_XII= sortrows(table(Vars',round(cf',2),round((pv'),3)),2,'descend','ComparisonMethod','abs');

% Generating Fig 18
f1=figure();
scatter(Ttct,Tfrac,'*black');
box on;
hold on;
yline(1,'-blue','LineWidth',1.5);
xlabel('$T_{s}^{5T2C5T}$','Interpreter','latex');
ylabel('$\frac{T_{s}^{1C5T1C5T}}{T_{s}^{5T2C5T}}$','Interpreter','latex');
xlim([129 12000]);
ylim([0.6 2.05]);
xticks([130 6000 12000]);
yticks([0.6 1 1.5 2]);
fontsize(f1,14,'point');
fontname(f1,"Arial");

% Generating Fig 19 a
f3=figure();
scatter(r1val,Tfrac,'*black');
box on;
hold on;
yline(1,'-blue','LineWidth',1.5);
xlabel('$r_{1}$','Interpreter','latex');
ylabel('$\frac{T_{s}^{1C5T1C5t}}{T_{s}^{5T2C5T}}$','Interpreter','latex');
xlim([0.001 0.025]);
xticks([0.001 0.012 0.025]);
yticks([0.6 1 2]);
fontsize(f3,14,'point');
fontname(f3,"Arial");

% Generating Fig 19 b
f4=figure();
scatter(rho4val,Tfrac,'*black');
box on;
hold on;
yline(1,'-blue','LineWidth',1.5);
xlabel('$\rho_{4}$','Interpreter','latex');
ylabel('$\frac{T_{s}^{1C5T1C5t}}{T_{s}^{5T2C5T}}$','Interpreter','latex');
xlim([0.01 0.1]);
xticks([0.01 0.05 0.1]);
ylim([0.6 2.1]);
yticks([0.6 1 2]);
fontsize(f4,14,'point');
fontname(f4,"Arial");
