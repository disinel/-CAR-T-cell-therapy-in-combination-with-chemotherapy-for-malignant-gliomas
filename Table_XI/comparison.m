clear variables;

load('CAR-T_TMZ_CAR-T_surv_max_r1_0_025_rho4_10_v_5_N=10000');

Tctc=Tsurv;

clearvars -except Tctc;

load('CAR-T_TMZ_CAR-T_TMZ_surv_max_r1_0_025_rho4_10_v_5_N=10000');

Tctct=Tsurv;

clearvars -except Tctc Tctct;


load('CAR-T_TMZ_Tsurv_max_r1_0_025_rho4_10_N=10000');

Tct=Tsurv;

clearvars -except Tctc Tctct Tct;

load('TMZ_CAR-T_TMZ_CAR-T_Tsurv_max_r1_0_025_rho4_10_N=10000');

Ttctc=Tsurv;

clearvars -except Tctc Tctct Tct Ttctc;


load('TMZ_CAR-T_TMZ_Tsurv_max_r1_0_025_rho4_10_v_5_N=10000');

Ttct=Tsurv;

clearvars -except Tctc Tctct Tct Ttctc Ttct Ttct;


load('TMZ_CAR-T_Tsurv_max_r1_0_025_rho4_10_N=10000');

Ttc=Tsurv;

clearvars -except Tctc Tctct Tct Ttctc Ttct Ttct Ttc;

Ts=[Ttct;Tct;Tctct;Ttctc;Ttc;Tctc];

for i=1:6
  for j=1:6
         [r,p]=corrcoef(Ts(i,:),Ts(j,:));
         R(i,j)=r(1,2); P(i,j)=p(1,2);
  end
end

M=zeros(2,length(Tctc));

for i=1:length(Tctc)
    [a,b]=max(Ts(:,i));
    M(1,i)=a; M(2,i)=b;
end


K={};

for i=1:length(Tctc)
    [a,b]=max(Ts(:,i));
    K(1,i)={a}; str=num2str(b);
    %arr=b;
    for j=1:6
        if j~=b
            if 100*(abs(a-Ts(j,i))/a)<10
                str=append(str,' ',num2str(j));
                %arr=[arr,j];
            end
        end
    end
    %K(2,i)=num2cell(arr);
     K(2,i)=cellstr(str);

end

Indx=zeros(1,length(Tct));

for i=1:length(Tct)
    arr=str2num(K{2,i});
    if length(arr)<=3
        Indx(i)=i;
    end
end

indx=find(Indx>0);

l6=0; l5=0; l4=0; l3=0; l2=0; l1=0;
for i=1:length(Tct)
    arr=str2num(K{2,i});
    if length(arr)==6
        l6=l6+1;
    end
    if length(arr)==5
        l5=l5+1;
    end
    if length(arr)==4
        l4=l4+1;
    end
    if length(arr)==3
        l3=l3+1;
    end
    if length(arr)==2
        l2=l2+1;
    end
    if length(arr)==1
        l1=l1+1;
    end
end

suml=l1+l2+l3+l4+l5+l6;

load('parameters_N=10000');

r1val_1=r1val(indx);
r1m=median(r1val);
r1m_1=median(r1val_1);
r1sh=-100*(r1m-r1m_1)/r1m;

alpha1val_1=alpha1val(indx);
alpha1m=median(alpha1val);
alpha1m_1=median(alpha1val_1);
alpha1sh=-100*(alpha1m-alpha1m_1)/alpha1m;

rho1val_1=rho1val(indx);
rho1m=median(rho1val);
rho1m_1=median(rho1val_1);
rho1sh=-100*(rho1m-rho1m_1)/rho1m;

rho2val_1=rho2val(indx);
rho2m=median(rho2val);
rho2m_1=median(rho2val_1);
rho2sh=-100*(rho2m-rho2m_1)/rho2m;

rho3val_1=rho3val(indx);
rho3m=median(rho3val);
rho3m_1=median(rho3val_1);
rho3sh=-100*(rho3m-rho3m_1)/rho3m;

epsilon1val_1=epsilon1val(indx);
epsilon1m=median(epsilon1val);
epsilon1m_1=median(epsilon1val_1);
epsilon1sh=-100*(epsilon1m-epsilon1m_1)/epsilon1m;

delta1val_1=delta1val(indx);
delta1m=median(delta1val);
delta1m_1=median(delta1val_1);
delta1sh=-100*(delta1m-delta1m_1)/delta1m;


T0val_1=T0val(indx);
T0m=median(T0val);
T0m_1=median(T0val_1);
T0sh=-100*(T0m-T0m_1)/T0m;

rho4val_1=rho4val(indx);
rho4m=median(rho4val);
rho4m_1=median(rho4val_1);
rho4sh=-100*(rho4m-rho4m_1)/rho4m;

delta2val_1=delta2val(indx);
delta2m=median(delta2val);
delta2m_1=median(delta2val_1);
delta2sh=-100*(delta2m-delta2m_1)/delta2m;


Vars={'r1','alpha1','rho1','rho2','rho3','epsilon1','delta1','T0','rho4','delta2'};

Shifts=[r1sh, alpha1sh, rho1sh, rho2sh, rho3sh, epsilon1sh, delta1sh, T0sh, rho4sh, delta2sh];

Table_XI= sortrows(table(Vars',round(Shifts,0)'),2,'descend','ComparisonMethod','abs');

disp(Table_XI)
