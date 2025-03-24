clear variables;
close all;


% Loading data from the file 'Notreatment_r2=2r1_max_rho4_10_N=10000' 
load('Notreatment_r2=2r1_max_rho4_10_N=10000');

Tsnt=Tsurv;

clearvars -except Tsnt;

% Loading data from the file 'CAR-T_trial_r2=2r1_max_rho4_10_v_5_N=10000_'
load('CAR-T_trial_r2=2r1_max_rho4_10_v_5_N=10000_')

% Computation of the correlation coefficient presented in Table XVI.  
[R,P]=corrcoef(r1val,Tsurv); cf(1)=R(1,2); pv(1)=P(1,2);   [R,P]=corrcoef(alpha1val,Tsurv); cf(2)=R(1,2); pv(2)=P(1,2);  

[R,P]=corrcoef(alpha3val,Tsurv); cf(3)=R(1,2);  pv(3)=P(1,2); [R,P]=corrcoef(epsilon1val,Tsurv); cf(4)=R(1,2);  pv(4)=P(1,2);

[R,P]=corrcoef(rho1val,Tsurv); cf(5)=R(1,2); pv(5)=P(1,2);  [R,P]=corrcoef(rho2val,Tsurv); cf(6)=R(1,2);  pv(6)=P(1,2);

[R,P]=corrcoef(rho3val,Tsurv);  cf(7)=R(1,2); pv(7)=P(1,2); [R,P]=corrcoef(rho4val,Tsurv); cf(8)=R(1,2); pv(8)=P(1,2);

[R,P]=corrcoef(T0val,Tsurv);  cf(9)=R(1,2); pv(9)=P(1,2); [R,P]=corrcoef(delta1val,Tsurv); cf(10)=R(1,2); pv(10)=P(1,2); 

[R,P]=corrcoef(delta2val,Tsurv); cf(11)=R(1,2);  pv(11)=P(1,2);


Vars={'r1','alpha1','alpha3','epsilon1','rho1','rho2','rho3','rho4','T0','delta1','delta2'};

Table_XVI= sortrows(table(Vars',round(cf',2),round((pv'),3)),2,'descend','ComparisonMethod','abs');


% Tsd represents the improvement in survival time for 10TMZ
% compared to no treatment
Tsd=Tsurv./Tsnt;

% Correlation coefficients for Table XVII
[R,P]=corrcoef(r1val,Tsd); cf(1)=R(1,2); pv(1)=P(1,2);   [R,P]=corrcoef(alpha1val,Tsd); cf(2)=R(1,2); pv(2)=P(1,2);  

[R,P]=corrcoef(alpha3val,Tsd); cf(3)=R(1,2);  pv(3)=P(1,2); [R,P]=corrcoef(epsilon1val,Tsd); cf(4)=R(1,2);  pv(4)=P(1,2);

[R,P]=corrcoef(rho1val,Tsd); cf(5)=R(1,2); pv(5)=P(1,2);  [R,P]=corrcoef(rho2val,Tsd); cf(6)=R(1,2);  pv(6)=P(1,2);

[R,P]=corrcoef(rho3val,Tsd);  cf(7)=R(1,2); pv(7)=P(1,2); [R,P]=corrcoef(rho4val,Tsd); cf(8)=R(1,2); pv(8)=P(1,2);

[R,P]=corrcoef(T0val,Tsd);  cf(9)=R(1,2); pv(9)=P(1,2); [R,P]=corrcoef(delta1val,Tsd); cf(10)=R(1,2); pv(10)=P(1,2); 

[R,P]=corrcoef(delta2val,Tsd); cf(11)=R(1,2);  pv(11)=P(1,2);


Vars={'r1','alpha1','alpha3','epsilon1','rho1','rho2','rho3','rho4','T0','delta1','delta2'};

Table_XVII= sortrows(table(Vars',round(cf',2),round((pv'),3)),2,'descend','ComparisonMethod','abs');
