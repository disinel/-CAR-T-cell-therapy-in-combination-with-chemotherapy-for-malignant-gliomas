clear variables;
close all;

% Loading parameters corresponding to a virtual population of 10,000 patients.
load('CAT-T_only_v_dep_N=10000');


f=figure();
plot(vv,Tsm,'-oblack','LineWidth',1.5);
xlabel('2v');
ylabel('$\widetilde{T}_{s}^{2C}$','interpreter','latex');
xticks([2*10^(7) 10^(9) 2*10^(9)]);
yticks([280  310 340]);
fontsize(f,14,'point');
fontname(f,"Arial");

