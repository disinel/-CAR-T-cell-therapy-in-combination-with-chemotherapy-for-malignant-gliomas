clear variables;
close all;

% Loading parameters corresponding to a virtual population of 10,000 patients.
load('CAT-T_only_L2_distributed_dep_N=10000');



f=figure();
plot(L2q,Tsm,'-oblack','LineWidth',1.5);
xlabel('L_{2}');
ylabel('$\widetilde{T}_{s}^{L_{2}C}$','interpreter','latex');
xticks([1 5 10]);
yticks([302  310 318]);
xlim([1 10]);
fontsize(f,14,'point');
fontname(f,"Arial");




