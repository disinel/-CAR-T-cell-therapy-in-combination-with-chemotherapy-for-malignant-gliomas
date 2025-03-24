clear variables;
close all;

% Load data from computations
load('CAT-T_only_T2_dep_N=10000');


f=figure();
plot(TT2,Tsm,'-oblack','LineWidth',1.5);
xlabel('T_{gap}^{(2)}');
ylabel('$\widetilde{T}_{s}^{2C}$','interpreter','latex');
xlim([1 30]);
ylim([308  316]);
xticks([1 15 30]);
yticks([308  312 316]);
fontsize(f,14,'point');
fontname(f,"Arial");

