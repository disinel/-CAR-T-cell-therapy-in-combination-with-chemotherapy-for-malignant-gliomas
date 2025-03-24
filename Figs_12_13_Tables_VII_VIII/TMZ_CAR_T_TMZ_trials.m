clear variables;
close all;

% Fixed parameters
g1=10^(10); g2=10^(10); g3=2*10^(9); v=1*10^(9); K=5*10^(12); mu=8.32; E0=1;

V=v/K; alpha2=2.5*10^(-10); a2=alpha2*K;  gamma1=g1/K; gamma2=g2/K; gamma3=g3/K;

% End time of the integration
Tfinal=20000; 

% The time interval between the combined application of TMZ and CAR-T cells (Tgap1), 
% and the time interval between two consecutive CAR-T cell applications (Tgap2).
Tgap1=1; Tgap2=7;

% Length and number of TMZ cycles
T1=28; L1=5; T11=1; L11=4;  LL1=5;

% Length and number of CAR-T cycles
T2=7; L2=2;

% Number of stored time points during ODEs integration
Np=10;

% Loading parameters corresponding to a virtual population of 10,000 patients.
load('parameters_N=10000');

% Array of survival times for each virtual patient.  
Tsurv=zeros(1,N);

tic

% Main loop computing survival time for each patient
parfor n=1:N 
 
Tc0=T0val(n); fs2=delta2val(n); frc=delta1val(n);  S10=(1-fs2-frc)*Tc0; S20=fs2*Tc0; RC0=frc*Tc0; 

r1=r1val(n); r2=0.5*r1; alpha1=alpha1val(n);  alpha3=alpha3val(n); rho1=rho1val(n); rho2=rho2val(n); rho3=rho3val(n); rho4=rho4val(n); epsilon1=epsilon1val(n);

[te] = calc_1cycle(r1,r2,alpha1,a2,alpha3,rho1,rho2,rho3,rho4,epsilon1,gamma1,gamma2,gamma3,mu,V,K,Np,T1,T11,L1,LL1,L11,T2,L2,S10,S20,RC0,E0,Tfinal,Tgap1,Tgap2);

Tsurv(n)=te;

end

toc


% Median survival time
median(Tsurv)


% Computation of the correlation coefficients
[R,P]=corrcoef(r1val,Tsurv); cf(1)=R(1,2); pv(1)=P(1,2);   [R,P]=corrcoef(alpha1val,Tsurv); cf(2)=R(1,2); pv(2)=P(1,2);  

[R,P]=corrcoef(alpha3val,Tsurv); cf(3)=R(1,2);  pv(3)=P(1,2); [R,P]=corrcoef(epsilon1val,Tsurv); cf(4)=R(1,2);  pv(4)=P(1,2);

[R,P]=corrcoef(rho1val,Tsurv); cf(5)=R(1,2); pv(5)=P(1,2);  [R,P]=corrcoef(rho2val,Tsurv); cf(6)=R(1,2);  pv(6)=P(1,2);

[R,P]=corrcoef(rho3val,Tsurv);  cf(7)=R(1,2); pv(7)=P(1,2); [R,P]=corrcoef(rho4val,Tsurv); cf(8)=R(1,2); pv(8)=P(1,2);

[R,P]=corrcoef(T0val,Tsurv);  cf(9)=R(1,2); pv(9)=P(1,2); [R,P]=corrcoef(delta1val,Tsurv); cf(10)=R(1,2); pv(10)=P(1,2); 

[R,P]=corrcoef(delta2val,Tsurv); cf(11)=R(1,2);  pv(11)=P(1,2);


Vars={'r1','alpha1','alpha3','epsilon1','rho1','rho2','rho3','rho4','T0','delta1','delta2'};

Table_VII= sortrows(table(Vars',round(cf',2),round((pv'),3)),2,'descend','ComparisonMethod','abs');





%  Saving generated data
save("TMZ_CAR-T_TMZ_Tsurv_max_r1_0_025_rho4_"+ceil(100*max(rho4val))+"_v_"+v/10^(8)+"_N="+N+".mat");



% This function computes the solution of the governing system of ODEs for the 
% given initial conditions. 
% If the total tumor size reaches 10^{12} (the number of cancer cells  
% considered fatal for the patient), the function stops computations 
% and returns the time moment when 10^{12} is reached, which is the 
% survival time. Otherwise, the corresponding variable 'te' is empty.

function [ t, y, te, ye] = calc_event(ics,T1,T2,r1,r2,alpha1,a2,alpha3,rho1,rho2,rho3,rho4,epsilon1,gamma1,gamma2,gamma3,mu,Np,K)
tspan = linspace(T1,T2,Np);
yy0=ics;
options = odeset('RelTol',1e-10,'Events',@events,NonNegative=[1,2,3,4,5]);

    function [position,isterminal,direction] = events(~,y)
        position = K.*(y(1)+y(2)+y(3))-10^(12);                                           
        isterminal = 1;                                                        
        direction = 0;                                                         
    end



[t,y, te,ye] = ode78(@func,tspan,yy0,options);




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

function [te] = calc_1cycle(r1,r2,alpha1,a2,alpha3,rho1,rho2,rho3,rho4,epsilon1,gamma1,gamma2,gamma3,mu,V,K,Np,T1,T11,L1,LL1,L11,T2,L2,S10,S20,RC0,E0,Tfinal,Tgap1,Tgap2)

ics(1)=S10/K; ics(2)=S20/K; ics(3)=RC0/K; ics(4)=0*10^(9)/K; ics(5)=E0;

Sy=(L1*(L11+1))*Np+Np+(L2-1)*Np+Np+(LL1*(L11+1))*Np+Np; 
 
Y=zeros(Sy,5); tt=zeros(Sy,1);



M=0;

for i=1:L1

 for j=1:L11

        [t,y,te]=calc_event(ics,T11*(j-1)+T1*(i-1),T11*j+T1*(i-1),r1,r2,alpha1,a2,alpha3,rho1,rho2,rho3,rho4,epsilon1,gamma1,gamma2,gamma3,mu,Np,K);
        if not(isempty(te)) return; end
        ics(1)=y(Np,1); ics(2)=y(Np,2); ics(3)=y(Np,3); ics(4)=y(Np,4); ics(5)=y(Np,5)+E0;
        tt(M*Np+1:M*Np+Np)=t;
        Y(M*Np+1:M*Np+Np,:)=y;
        M=M+1;
end


[t,y,te]=calc_event(ics,t(length(t)),T1*i,r1,r2,alpha1,a2,alpha3,rho1,rho2,rho3,rho4,epsilon1,gamma1,gamma2,gamma3,mu,Np,K);
if not(isempty(te)) return; end
tt(M*Np+1:M*Np+Np)=t;
Y(M*Np+1:M*Np+Np,:)=y;
M=M+1;
   
ics(1)=y(Np,1); ics(2)=y(Np,2); ics(3)=y(Np,3); ics(4)=y(Np,4); ics(5)=y(Np,5)+E0;

end



ics(1)=Y(M*Np,1); ics(2)=Y(M*Np,2); ics(3)=Y(M*Np,3); ics(4)=Y(M*Np,4); ics(5)=Y(M*Np,5);

[t,y,te]=calc_event(ics,L1*T1,L1*T1+Tgap1,r1,r2,alpha1,a2,alpha3,rho1,rho2,rho3,rho4,epsilon1,gamma1,gamma2,gamma3,mu,Np,K);
if not(isempty(te)) return; end
tt(M*Np+1:M*Np+Np)=t;
Y(M*Np+1:M*Np+Np,:)=y;
M=M+1;

ics(1)=Y(M*Np,1); ics(2)=Y(M*Np,2); ics(3)=Y(M*Np,3); ics(4)=Y(M*Np,4)+V; ics(5)=Y(M*Np,5);

for i=1:L2-1
    [t,y,te]=calc_event(ics,L1*T1+Tgap1+T2*(i-1),L1*T1+Tgap1+T2*i,r1,r2,alpha1,a2,alpha3,rho1,rho2,rho3,rho4,epsilon1,gamma1,gamma2,gamma3,mu,Np,K);
    if not(isempty(te)) return; end
    ics(1)=y(Np,1); ics(2)=y(Np,2); ics(3)=y(Np,3); ics(4)=y(Np,4)+V; ics(5)=y(Np,5);
    tt(M*Np+1:M*Np+Np)=t;
    Y(M*Np+1:M*Np+Np,:)=y;
    M=M+1;
end



[t,y,te]=calc_event(ics,L1*T1+Tgap1+(L2-1)*T2,L1*T1+Tgap1+(L2-1)*T2+Tgap2,r1,r2,alpha1,a2,alpha3,rho1,rho2,rho3,rho4,epsilon1,gamma1,gamma2,gamma3,mu,Np,K);
if not(isempty(te)) return; end
tt(M*Np+1:M*Np+Np)=t;
Y(M*Np+1:M*Np+Np,:)=y;
M=M+1;

ics(1)=Y(M*Np,1); ics(2)=Y(M*Np,2); ics(3)=Y(M*Np,3); ics(4)=Y(M*Np,4); ics(5)=Y(M*Np,5)+E0;

for i=1:LL1

 for j=1:L11

     [t,y,te]=calc_event(ics,L1*T1+Tgap1+(L2-1)*T2+Tgap2+T11*(j-1)+T1*(i-1),L1*T1+Tgap1+(L2-1)*T2+Tgap2+T11*j+T1*(i-1),  r1,r2,alpha1,a2,alpha3,rho1,rho2,rho3,rho4,epsilon1,gamma1,gamma2,gamma3,mu,Np,K);
     if not(isempty(te)) return; end
     ics(1)=y(Np,1); ics(2)=y(Np,2); ics(3)=y(Np,3); ics(4)=y(Np,4); ics(5)=y(Np,5)+E0;
     tt(M*Np+1:M*Np+Np)=t;
     Y(M*Np+1:M*Np+Np,:)=y;
     M=M+1;
end


[t,y,te]=calc_event(ics,t(length(t)),L1*T1+Tgap1+(L2-1)*T2+Tgap2+T1*i,r1,r2,alpha1,a2,alpha3,rho1,rho2,rho3,rho4,epsilon1,gamma1,gamma2,gamma3,mu,Np,K);
if not(isempty(te)) return; end
tt(M*Np+1:M*Np+Np)=t;
Y(M*Np+1:M*Np+Np,:)=y;
M=M+1;
   
ics(1)=y(Np,1); ics(2)=y(Np,2); ics(3)=y(Np,3); ics(4)=y(Np,4); ics(5)=y(Np,5)+E0;

end


ics(1)=Y(M*Np,1); ics(2)=Y(M*Np,2); ics(3)=Y(M*Np,3); ics(4)=Y(M*Np,4); ics(5)=Y(M*Np,5);

[t,y,te]= calc_event(ics,L1*T1+Tgap1+(L2-1)*T2+Tgap2+T1*LL1,Tfinal,r1,r2,alpha1,a2,alpha3,rho1,rho2,rho3,rho4,epsilon1,gamma1,gamma2,gamma3,mu,Np,K);
tt(M*Np+1:M*Np+length(t))=t;
Y(M*Np+1:M*Np+length(t),:)=y;




end
