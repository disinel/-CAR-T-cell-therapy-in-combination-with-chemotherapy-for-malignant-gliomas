clear variables;
close all;

% Loading data from the file 'CTCT_n=1454_'
load('CTCT_n=1454_');

% Generating Fig. 17, bottom panel a
f=figure();
plot(tt,YY(:,1),'-','color','#77AC30','LineWidth',1.5);
hold on;
plot(tt,YY(:,2),'-','color','red','LineWidth',1.5);
hold on;
plot(tt,YY(:,3),'-','color','#EDB120','LineWidth',1.5);
hold on;
plot(tt,Tc,'-black','LineWidth',1.5);
hold on;
xticks([0 150 300]);
xlabel('t');
fontsize(f,14,'point');
fontname(f,"Arial");

clearvars -except f;

% Loading data from the file 'TCT_n=1454_'
load('TCT_n=1454_');

% Generating Fig. 17, bottom panel a
plot(tt,YY(:,1),'--','color','#77AC30','LineWidth',1.5);
hold on;
plot(tt,YY(:,2),'--','color','red','LineWidth',1.5);
hold on;
plot(tt,YY(:,3),'--','color','#EDB120','LineWidth',1.5);
hold on;
plot(tt,Tc,'--black','LineWidth',1.5);
xticks([0 125 250]);
yticks([1*10^(10) 1.25*10^(11) 2.5*10^(11)]);
xlim([0 250]);
xlabel('t');
legend('S','R_{C}','R_{E}','T','S','R_{C}','R_{E}','T','location','northwest');
fontsize(f,14,'point');
fontname(f,"Arial");
