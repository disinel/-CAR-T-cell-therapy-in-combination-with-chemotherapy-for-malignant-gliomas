clear variables;
close all;


% Loading parameters 
load('Notreatment_max_rho4_10_N=10000');

Tsnt=Tsurv;

clearvars -except Tsnt;

load('CAR-T_trial_max_rho4_10_v_5_N=10000_v=0_5_10_9');

% Computation of the correlation coefficient presented in Table IV.  

[R,P]=corrcoef(r1val,Tsurv); cf(1)=R(1,2); pv(1)=P(1,2);   [R,P]=corrcoef(rho1val,Tsurv); cf(2)=R(1,2); pv(2)=P(1,2);  [R,P]=corrcoef(rho2val,Tsurv); cf(3)=R(1,2);  pv(3)=P(1,2);

[R,P]=corrcoef(rho3val,Tsurv);  cf(4)=R(1,2); pv(4)=P(1,2); [R,P]=corrcoef(rho4val,Tsurv); cf(5)=R(1,2); pv(5)=P(1,2); [R,P]=corrcoef(T0val,Tsurv); cf(6)=R(1,2);  pv(6)=P(1,2);

[R,P]=corrcoef(delta1val,Tsurv); cf(7)=R(1,2); pv(7)=P(1,2);  [R,P]=corrcoef(delta2val,Tsurv); cf(8)=R(1,2); pv(8)=P(1,2);

Vars={'r1','rho1','rho2','rho3','rho4','T0','delta1','delta2'};

Table_IV= sortrows(table(Vars',round(cf',2),round((pv'),3)),2,'descend','ComparisonMethod','abs');


% Tsd represents the improvement in survival time
% for 10TMZ compared to no treatment.  

Tsd=Tsurv./Tsnt;

% Correlation coefficients for Table V
[R,P]=corrcoef(r1val,Tsd); cf(1)=R(1,2); pv(1)=P(1,2);   [R,P]=corrcoef(alpha1val,Tsd); cf(2)=R(1,2); pv(2)=P(1,2);  

[R,P]=corrcoef(alpha3val,Tsd); cf(3)=R(1,2);  pv(3)=P(1,2); [R,P]=corrcoef(epsilon1val,Tsd); cf(4)=R(1,2);  pv(4)=P(1,2);

[R,P]=corrcoef(rho1val,Tsd); cf(5)=R(1,2); pv(5)=P(1,2);  [R,P]=corrcoef(rho2val,Tsd); cf(6)=R(1,2);  pv(6)=P(1,2);

[R,P]=corrcoef(rho3val,Tsd);  cf(7)=R(1,2); pv(7)=P(1,2); [R,P]=corrcoef(rho4val,Tsd); cf(8)=R(1,2); pv(8)=P(1,2);

[R,P]=corrcoef(T0val,Tsd);  cf(9)=R(1,2); pv(9)=P(1,2); [R,P]=corrcoef(delta1val,Tsd); cf(10)=R(1,2); pv(10)=P(1,2); 

[R,P]=corrcoef(delta2val,Tsd); cf(11)=R(1,2);  pv(11)=P(1,2);


Vars={'r1','alpha1','alpha3','epsilon1','rho1','rho2','rho3','rho4','T0','delta1','delta2'};

Table_V= sortrows(table(Vars',round(cf',2),round((pv'),3)),2,'descend','ComparisonMethod','abs');


% Generating Fig. 11
f=figure();
scatter(rho4val,Tsd,'*black');
box on;
hold on;
xlabel('$\rho_{4}$','Interpreter','latex','FontSize',14);
ylabel('$\frac{T_{s}^{2C}}{T_{s}^{NT}}$','Interpreter','latex');
axis([0.002 0.1 1 1.9]);
xticks([0.002 0.05 0.1]);
yticks([1 1.5 1.9]);
fontsize(f,14,"points");
fontname(f,'arial');

