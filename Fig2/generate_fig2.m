close all;
clear variables;

% Load data from computations

load('TMZ_only_Tsurv_r2=0_5_r1_max_r1_0_025_N=10000');

for q=1:Q
    Tsm1(q)=median(Tsurv(q,:));
end

load('TMZ_only_Tsurv_r2=r1_max_r1_0_025_N=10000');

for q=1:Q
    Tsm2(q)=median(Tsurv(q,:));
end

load('TMZ_only_Tsurv_r2=2r1_max_r1_0_025_N=10000');

for q=1:Q
    Tsm3(q)=median(Tsurv(q,:));
end

% Addition of median survival times for the absence of treatment and
% for r_{2}=r_{1}/2, r_{2}=r_{1} and r_{2}=2r_{1}, respectively
Tsm1=[268,Tsm1];
Tsm2=[264,Tsm2];
Tsm3=[221,Tsm3];
L1d=[0,L1d];


f=figure();
plot(L1d,Tsm1,'o-black','LineWidth',1.5);
hold on;
plot(L1d,Tsm2,'o--blue','LineWidth',1.5);
hold on;
plot(L1d,Tsm3,'o-.black','LineWidth',1.5);
legend('r_{2}=0.5r_{1}','r_{2}=r_{1}','r_{2}=2r_{1}','location','northwest');
xlabel('L_{1}');
ylabel('$\widetilde{T}_{s}^{L_{1}T}$','Interpreter','latex');
xlim([0 10]);
ylim([160 580]);
xticks([0 5 10]);
yticks([160 300 440 580]);
fontsize(f,14,'points');





