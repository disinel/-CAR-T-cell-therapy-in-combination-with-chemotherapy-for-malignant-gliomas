clear variables;
close all;

% Fixed parameters
g1=10^(10); g2=10^(10); g3=2*10^(9); v=0.5*10^(9); K=5*10^(12); mu=8.32; E0=1;

V=v/K; alpha2=2.5*10^(-10); a2=alpha2*K;  gamma1=g1/K; gamma2=g2/K; gamma3=g3/K;

% End time of the integration
Tfinal=20000; 



% Length and number of CAR-T cycles and
% Number of stored time points during ODEs integration
T2=7; L2=2; Np=10;


tic
% Loading parameters corresponding to a virtual population of 10,000 patients.
load('parameters_N=10000');

% Array of survival times for each virtual patient.  
Tsurv=zeros(1,N);

% Main loop computing survival time for each patient
parfor n=1:N

      Tc0=T0val(n); fs2=delta2val(n); frc=delta1val(n);  S10=(1-fs2-frc)*Tc0; S20=fs2*Tc0; RC0=frc*Tc0; 

      r1=r1val(n); r2=2*r1; alpha1=alpha1val(n);  alpha3=alpha3val(n); rho1=rho1val(n); rho2=rho2val(n); rho3=rho3val(n); rho4=rho4val(n); epsilon1=epsilon1val(n);
    
      [te] = calc_1cycle(r1,r2,alpha1,a2,alpha3,rho1,rho2,rho3,rho4,epsilon1,gamma1,gamma2,gamma3,mu,V,K,Np,T2,L2,S10,S20,RC0,E0,Tfinal);

      Tsurv(n)=te;

 end

toc

% Computation of the correlation coefficient presented in Table XVI and XVII.  

[R,P]=corrcoef(r1val,Tsurv); cf(1)=R(1,2); pv(1)=P(1,2);   [R,P]=corrcoef(rho1val,Tsurv); cf(2)=R(1,2); pv(2)=P(1,2);  [R,P]=corrcoef(rho2val,Tsurv); cf(3)=R(1,2);  pv(3)=P(1,2);

[R,P]=corrcoef(rho3val,Tsurv);  cf(4)=R(1,2); pv(4)=P(1,2); [R,P]=corrcoef(rho4val,Tsurv); cf(5)=R(1,2); pv(5)=P(1,2); [R,P]=corrcoef(T0val,Tsurv); cf(6)=R(1,2);  pv(6)=P(1,2);

[R,P]=corrcoef(delta1val,Tsurv); cf(7)=R(1,2); pv(7)=P(1,2);  [R,P]=corrcoef(delta2val,Tsurv); cf(8)=R(1,2); pv(8)=P(1,2);

Vars={'r1','rho1','rho2','rho3','rho4','T0','delta1','delta2'};

Table = table(Vars',round(abs(cf'),2),round(abs(pv'),5));

Table1 = sortrows(Table,2,'descend');

% Median survival time
median(Tsurv)


Tm=max(Tsurv); Tmct=ceil(Tm);

for i=1:Tmct
    Nsurv_ct(i)=length(Tsurv(Tsurv>=i));
end

Tl_ct=linspace(1,Tmct,Tmct);

%  Saving generated data
save("CAR-T_trial_r2=2r1_max_rho4_"+ceil(100*max(rho4val))+"_v_"+v/10^(8)+"_N="+N+"_");




% This function computes the solution of the governing system of ODEs for the 
% given initial conditions. 
% If the total tumor size reaches 10^{12} (the number of cancer cells  
% considered fatal for the patient), the function stops computations 
% and returns the time moment when 10^{12} is reached, which is the 
% survival time. Otherwise, the corresponding variable 'te' is empty.



function [ t, y, te] = calc_event(ics,T1,T2,r1,r2,alpha1,a2,alpha3,rho1,rho2,rho3,rho4,epsilon1,gamma1,gamma2,gamma3,mu,Np,K)
tspan = linspace(T1,T2,Np);
yy0=ics;
options = odeset('RelTol',1e-10,'Events',@events,NonNegative=[1,2,3,4,5]);

    function [position,isterminal,direction] = events(~,y)
        position = K.*(y(1)+y(2)+y(3))-10^(12);                                           
        isterminal = 1;                                                        
        direction = 0;                                                         
    end



[t,y, te] = ode78(@func,tspan,yy0,options);




function dy = func(~,y)
      dy=zeros(5,1);
      dy(1)=(r1.*(1-y(1)-y(2)-y(3))-(alpha1+epsilon1)*y(5)-a2*y(4))*y(1);
      dy(2)=(r1.*(1-y(1)-y(2)-y(3))-(alpha1+epsilon1)*y(5)).*y(2);
      dy(3)=r2*(1-y(1)-y(2)-y(3))*y(3)-a2*y(3)*y(4)+epsilon1*y(5)*(y(1)+y(2));
      dy(4)=-rho1*y(4)+(rho2.*y(1)*y(4))./(gamma1+y(1))+(rho3*y(3)*y(4))/(gamma2+y(3))-rho4*(y(1)+y(2)+y(3))*y(4)/(gamma3+y(4))-alpha3*y(5)*y(4);
      dy(5)=-mu*y(5);
     end

end

% This function, for a given virtual patient, follows the CAR-T cell application
% protocol and returns patient's survival time.

function [te] = calc_1cycle(r1,r2,alpha1,a2,alpha3,rho1,rho2,rho3,rho4,epsilon1,gamma1,gamma2,gamma3,mu,V,K,Np,T2,L2,S10,S20,RC0,E0,Tfinal)

ics(1)=S10/K; ics(2)=S20/K; ics(3)=RC0/K; ics(4)=V; ics(5)=0;

Sy=(L2-1)*Np+Np; 
 
Y=zeros(Sy,5); tt=zeros(Sy,1);

M=0;

for i=1:L2-1
    [t,y, te]= calc_event(ics,T2*(i-1),T2*i,r1,r2,alpha1,a2,alpha3,rho1,rho2,rho3,rho4,epsilon1,gamma1,gamma2,gamma3,mu,Np,K);
    if not(isempty(te)) return; end
    ics(1)=y(Np,1); ics(2)=y(Np,2); ics(3)=y(Np,3); ics(4)=y(Np,4)+V; ics(5)=y(Np,5);
    tt(M*Np+1:M*Np+Np)=t;
    Y(M*Np+1:M*Np+Np,:)=y;
    M=M+1;
end

[t,y,te]=calc_event(ics,(L2-1)*T2,(L2-1)*T2+Tfinal,r1,r2,alpha1,a2,alpha3,rho1,rho2,rho3,rho4,epsilon1,gamma1,gamma2,gamma3,mu,Np,K);
if not(isempty(te)) return; end
tt(M*Np+1:M*Np+Np)=t;
Y(M*Np+1:M*Np+Np,:)=y;
M=M+1;





end
