clear variables;
close all;

% Loading data from the file 'CTCT_n=3999_'
load('CTCT_n=3999_');

% Generating Fig. 17, upper panel b
f=figure();
colororder({'black','black'});
plot(tt,YY(:,4)./v,'-','color','red','LineWidth',1.5);
hold on;
plot(tt,Y(:,5),'-','color','#EDB120','LineWidth',1.5);
xlabel('t');
ylabel('$\frac{C}{v}$, E','interpreter','latex');

Tctct=Tc;
tt1=tt;


clearvars -except f tt1 Tctct;

% Loading data from the file 'TCT_n=3999_'
load('TCT_n=3999_');

% Generating Fig. 17, upper panel b
plot(tt,YY(:,4)./v,'--','color','red','LineWidth',1.5);
hold on;
plot(tt,Y(:,5),'--','color','#EDB120','LineWidth',1.5);
xticks([0 200 400]);
xlim([0 400]);
yticks([0 1.4 3]);
xlabel('t');
ylabel('$\frac{C}{v}$, E','interpreter','latex');

yyaxis right;
plot(tt1,Tctct,'-black','LineWidth',1.5);
hold on;
plot(tt,Tc,'--black','LineWidth',1.5);
yticks([1*10^(10) 1.25*10^(11) 2.5*10^(11)]);
legend('$\frac{C}{v}$','E','$\frac{C}{v}$','E','T','T','location','northwest','interpreter','latex');

fontsize(f,14,'point');
fontname(f,"Arial");


