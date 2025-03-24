clear variables;
close all;

% Fixed parameters
g1=10^(10); g2=10^(10); g3=2*10^(9); v=0.5*10^(9); K=5*10^(12); mu=8.32; E0=1;

V=v/K; alpha2=2.5*10^(-10); a2=alpha2*K;  gamma1=g1/K; gamma2=g2/K; gamma3=g3/K;

% End time of the integration
Tfinal=500;




tic

% Length of a TMZ cycle: the time interval in days between TMZ  
% applications. L11+1 represents the total number of TMZ applications per cycle. 
T1=28; L1=10; T11=1; L11=4; 

% Number of stored time points during ODEs integration and total number  
% of time points for the entire therapy.
Np=100; Sy=(L1*(L11+1))*Np+Np; 


% Loading parameters corresponding to a virtual population of 10,000 patients.
load('parameters_N=10000');

% Median values of the parameters
Tc0=median(T0val); fs2=median(delta2val); frc=median(delta1val);  S10=(1-fs2-frc)*Tc0; S20=fs2*Tc0; RC0=frc*Tc0; 

r1=median(r1val); r2=median(r2val); alpha1=median(alpha1val);  alpha3=median(alpha3val); rho1=median(rho1val); rho2=median(rho2val); rho3=median(rho3val); 
 
rho4=median(rho4val); epsilon1=median(epsilon1val);

% Computation of the solution of the ODE system corresponding to a given  
% protocol, in this case, 10 TMZ.
[tt, Y] = calc_1cycle(r1,r2,alpha1,a2,alpha3,rho1,rho2,rho3,rho4,epsilon1,gamma1,gamma2,gamma3,mu,V,K,Np,T1,T11,L1,L11,S10,S20,RC0,E0,Tfinal);

toc

% Conversion to dimensional variables
YY=K.*Y;

% Total tumor size
Tc=K.*(Y(:,1)+Y(:,2)+Y(:,3));






% Generating Fig. 4


f=figure();
plot(tt,YY(:,1),'-','color','#77AC30','LineWidth',1.5);
hold on;
plot(tt,YY(:,2),'-','color','red','LineWidth',1.5);
hold on;
plot(tt,YY(:,3),'-','color','#EDB120','LineWidth',1.5);
hold on;
plot(tt,Tc,'-black','LineWidth',1.5);
axis([min(tt) max(tt) min(YY(:,1)) max(Tc)]);
xticks([0 250 500]);
yticks([1*10^(10) 3.5*10^(11) 7*10^(11)]);
xlabel('t');
legend('S','R_{C}','R_{E}','T','location','northwest');
fontsize(f,14,'point');
fontname(f,"Arial");

f2=figure();
colororder({'black','black'});
plot(tt,Y(:,5),'-.','color','#EDB120','LineWidth',1.5);
xticks([0 250 500]);
yticks([0 1]);
ylim([0 1]);
xlabel('t');
ylabel('E','interpreter','latex');

yyaxis right;
plot(tt,Tc,'-black','LineWidth',1.5);
yticks([1*10^(10) 3.5*10^(11) 7*10^(11)]);
ylabel('T');
legend('E','T','location','northeast','interpreter','latex');

fontsize(f2,14,'point');
fontname(f2,"Arial");






% This function computes the solution of the governing system of ODEs  
% for given initial conditions over a specified time interval and for  
% specified parameter values.


function [ t, y] = calc(ics,T1,T2,r1,r2,alpha1,a2,alpha3,rho1,rho2,rho3,rho4,epsilon1,gamma1,gamma2,gamma3,mu,V,Np)
tspan = linspace(T1,T2,Np);
yy0=ics;
options = odeset('RelTol',1e-10,NonNegative=[1,2,3,4,5]); 
[t,y] = ode15s(@func,tspan,yy0,options);




function dy = func(~,y)
      dy=zeros(5,1);
      dy(1)=(r1.*(1-y(1)-y(2)-y(3))-(alpha1+epsilon1)*y(5)-a2*y(4))*y(1);
      dy(2)=(r1.*(1-y(1)-y(2)-y(3))-(alpha1+epsilon1)*y(5)).*y(2);
      dy(3)=r2*(1-y(1)-y(2)-y(3))*y(3)-a2*y(3)*y(4)+epsilon1*y(5)*(y(1)+y(2));
      dy(4)=0*V-rho1*y(4)+(rho2.*y(1)*y(4))./(gamma1+y(1))+(rho3*y(3)*y(4))/(gamma2+y(3))-rho4*(y(1)+y(2)+y(3))*y(4)/(gamma3+y(4))-alpha3*y(5)*y(4);
      dy(5)=-mu*y(5);
     end

end

% This function, for a given virtual patient, follows the TMZ application  
% protocol and returns the time points and the variables S, R_{C}, R_{E}, C and
% E

function [tt, Y] = calc_1cycle(r1,r2,alpha1,a2,alpha3,rho1,rho2,rho3,rho4,epsilon1,gamma1,gamma2,gamma3,mu,V,K,Np,T1,T11,L1,L11,S10,S20,RC0,E0,Tfinal)

ics(1)=S10/K; ics(2)=S20/K; ics(3)=RC0/K; ics(4)=0*10^(9)/K; ics(5)=E0;

Sy=(L1*(L11+1))*Np+Np; 
 
Y=zeros(Sy,5); tt=zeros(Sy,1);

M=0;

for i=1:L1

 for j=1:L11

        [t,y]=calc(ics,T11*(j-1)+T1*(i-1),T11*j+T1*(i-1),r1,r2,alpha1,a2,alpha3,rho1,rho2,rho3,rho4,epsilon1,gamma1,gamma2,gamma3,mu,V,Np);
        ics(1)=y(Np,1); ics(2)=y(Np,2); ics(3)=y(Np,3); ics(4)=y(Np,4); ics(5)=y(Np,5)+E0;
        tt(M*Np+1:M*Np+Np)=t;
        Y(M*Np+1:M*Np+Np,:)=y;
        M=M+1;
end


[t,y]=calc(ics,t(length(t)),T1*i,r1,r2,alpha1,a2,alpha3,rho1,rho2,rho3,rho4,epsilon1,gamma1,gamma2,gamma3,mu,V,Np);
tt(M*Np+1:M*Np+Np)=t;
Y(M*Np+1:M*Np+Np,:)=y;
M=M+1;
   
ics(1)=y(Np,1); ics(2)=y(Np,2); ics(3)=y(Np,3); ics(4)=y(Np,4); ics(5)=y(Np,5)+E0;

end



ics(1)=Y(M*Np,1); ics(2)=Y(M*Np,2); ics(3)=Y(M*Np,3); ics(4)=Y(M*Np,4); ics(5)=Y(M*Np,5);

[t,y]=calc(ics,L1*T1,Tfinal,r1,r2,alpha1,a2,alpha3,rho1,rho2,rho3,rho4,epsilon1,gamma1,gamma2,gamma3,mu,V,Np);
tt(M*Np+1:M*Np+Np)=t;
Y(M*Np+1:M*Np+Np,:)=y;





end
