%% This MATLAB code scans through the parameter sets for the HNSCC TME model in a resource-competition setting (Fig S1)
%%% Purpose of the study: The effects of varying resource competition on TME subtypes for varying maximum resource supply rate

%% Initialisation

P_para=2000*lhsdesign(24,10000);
y_0=5*ones(24,1); y_0(9)=0;
anti_PD1=2;
alpha1=[0.01 0.05 0.1 0.3];
tspan=[0 70000];
CAF_ST=[];
T_ST=[];
TUM_ST=[];
n=1;
X_ST=[];
X_ST_Post=[];
Tum_Prol_D_Bal=[];
T_Prol_D_Bal=[];
CAF_Prol_D_Bal=[];
Y_RM=10^3;                                             % Maximum resource amounts
for k=1:1:length(alpha1)
    P=HNSCC_parameters(alpha1(k));
    P(8)=0;
    P(11)=Y_RM;
    P(16)=1500;
    P(1:6)=P(1:6)/10;
for j=1:1:length(P_para(1,:))
    P(12)=0.15*P_para(1,j);
    P(14)=P_para(2,j);
    P(30)=P_para(3,j);
    P(32)=P_para(4,j);
    P(33)=P_para(5,j);
    P(36)=P_para(6,j);
    P(37)=P_para(7,j);
    P(38)=P_para(8,j);
    P(39)=P_para(9,j);
    P(40)=P_para(10,j);
    P(41)=P_para(11,j);
    P(42)=P_para(12,j);
    P(43)=P_para(13,j);
    P(44)=P_para(14,j);
    P(45)=P_para(15,j);
    P(48)=P_para(16,j);
    P(49)=P_para(17,j);
    P(50)=P_para(18,j);
    P(51)=P_para(19,j);
    P(53)=P_para(21,j);
    P(54)=P_para(22,j);
    P(55)=P_para(23,j);
    P(58)=P_para(24,j);
   [t_pre,x_pre]=ode23s(@(t,y)HNSCC_mod_Res_Comp(t,y,0,0,0,0,0,0,P),tspan,y_0);
   X_ST(k*j,:)=x_pre(length(x_pre),:);
Tum_Prol_D_Bal(k*j)=(1+P(12))*(1+P(14))*max(P(1:6))/(P(16)+P(26));
T_Prol_D_Bal(k*j)=P(30)*(1+P(37))*(1+P(38))/(P(36)+P(42));
CAF_Prol_D_Bal(k*j)=P(48)*(1+P(49))*(1+P(50))*(1+P(54))/(P(58)+P(56));
end 
for k1=1:1:1000
    P(12)=0.05*P_para(1,k1);
    P(30)=P_para(2,k1);
    P(32)=P_para(3,k1);
    P(33)=P_para(4,k1);
    P(36)=P_para(5,k1);
    P(37)=P_para(6,k1);
    P(38)=P_para(7,k1);
    P(39)=P_para(8,k1);
    P(40)=P_para(9,k1);
    P(41)=P_para(10,k1);
    P(42)=P_para(11,k1);
    P(43)=P_para(12,k1);
    P(44)=P_para(13,k1);
    P(45)=P_para(14,k1);
    P(58)=P_para(15,k1);
    P(47)=P(58)*lhsdesign(1,1);
    P(48)=10*lhsdesign(1,1);
    P(49)=10*lhsdesign(1,1);
    P(50)=10*lhsdesign(1,1);
    P(51)=10*lhsdesign(1,1);
    P(52)=P_para(19,k1);
    P(53)=P_para(20,k1);
    P(54)=10*lhsdesign(1,1);
    P(55)=10*lhsdesign(1,1);
   [t_pre,x_pre]=ode23s(@(t,y)HNSCC_mod_Res_Comp(t,y,0,0,0,0,0,0,P),tspan,y_0);
   X_ST(k*(j+k1),:)=x_pre(length(x_pre),:);
   Tum_Prol_D_Bal(k*(j+k1))=P(12)*P(14)*max(P(1:6))/(P(16)+P(26));
   T_Prol_D_Bal(k*(j+k1))=P(30)*P(37)/(P(40)+P(43))*P(38)/(P(36)+P(42));
   CAF_Prol_D_Bal(k*(j+k1))=P(48)*P(49)*P(50)*P(54)/(P(58)+P(56));
end 
end
Tum_PD=Tum_Prol_D_Bal;
TK_PD=T_Prol_D_Bal;
CAF_PD=CAF_Prol_D_Bal;
Bal_PD=[Tum_PD' TK_PD' CAF_PD'];

C_Tum=sum(X_ST(:,1:6),2)/30000;
C_Tk=sum(X_ST(:,8),2)/5000;
C_CAF=sum(X_ST(:,14),2)/5000;
ST=[C_Tum C_Tk C_CAF];
L1=find(ST(:,1)<0);
ST(L1,:)=[];
Bal_PD(L1,:)=[];
L2=find(ST(:,1)>1);
ST(L2,:)=[];
Bal_PD(L2,:)=[];
L3=find(ST(:,2)<0);
ST(L3,:)=[];
Bal_PD(L3,:)=[];
L4=find(ST(:,2)>1);
ST(L4,:)=[];
Bal_PD(L4,:)=[];
L5=find(ST(:,3)<0);
ST(L5,:)=[];
Bal_PD(L5,:)=[];
L6=find(ST(:,3)>1);
ST(L6,:)=[];
Bal_PD(L6,:)=[];

subplot(length(alpha1),length(Y_RM),k)

C=ST1;
%C=[sum(X_ST(1:6))/30000 X_ST(8)/5000 X_ST(14)/5000];
figure
ax3=axes();
S=20*ones(length(ST1),1);
scatter3(log(Bal_PD(:,1)),log(Bal_PD(:,2)),log(Bal_PD(:,3)),S,C,'filled')
xlabel('\zeta_{Tumor cell}','FontSize',20,'FontWeight','bold','FontName','Palatino Linotype')
ylabel('\zeta_{Killer T cell}','FontSize',20,'FontWeight','bold','FontName','Palatino Linotype')
zlabel('\zeta_{CAF}','FontSize',20,'FontWeight','bold','FontName','Palatino Linotype')
set(ax3,'FontName','Palatino Linotype','FontSize',19,'FontWeight','bold');

ST1=ST;
P_Count=find(ST(:,2)>0.1 & ST(:,1)>0.9 & ST(:,3)<0.5);
ST1(P_Count,:)=[];
C=ones(length(ST1),3);
FD_Count=find(ST1(:,1)<0.6 & ST1(:,2)>0.05 & ST1(:,3)<0.2);
IR_Count=find(ST1(:,1)<0.6 & ST1(:,2)>0.01 & ST1(:,3)>0.2);
FR_Count=find(ST1(:,1)>0.6 & ST1(:,2)>0.05 & ST1(:,3)>0.2);
ID_Count=find(ST1(:,1)>0.35 & ST1(:,2)<0.05 & ST1(:,3)>0.1);
 D_Count=find(ST1(:,1)>0.36 & ST1(:,2)<0.05 & ST1(:,3)<0.2);

for n_FD=1:1:length(FD_Count)
    C(FD_Count(n_FD),:)=[1 0.76 0.05];
end

for n_IR=1:1:length(IR_Count)
    C(IR_Count(n_IR),:)=[0.06 0.4 0.65];
end

for n_FR=1:1:length(FR_Count)
    C(FR_Count(n_FR),:)=[0.285 0.2 0.7];
end

for n_ID=1:1:length(ID_Count)
    C(ID_Count(n_ID),:)=[0.9 0.4 0];
end

for n_D=1:1:length(D_Count)
    C(D_Count(n_D),:)=[0.9 0.8 0.4];
end
ST1(ID_Count,:)=[];
C(ID_Count,:)=[];
figure
ax3=axes();
S=20*ones(length(IR_Count),1);
scatter3(ST1(IR_Count,1),ST1(IR_Count,2),ST1(IR_Count,3),S,C(IR_Count,:),'filled')
xlim([0 1])
ylim([0 1])
zlim([0 1])
xlabel('Tumor cell','FontSize',20,'FontWeight','bold','FontName','Palatino Linotype')
ylabel('Killer T cell','FontSize',20,'FontWeight','bold','FontName','Palatino Linotype')
zlabel('CAF','FontSize',20,'FontWeight','bold','FontName','Palatino Linotype')
set(ax3,'FontName','Palatino Linotype','FontSize',19,'FontWeight','bold');
%C=ST1;

figure
ax4=axes();
S1=20*ones(length(IR_Count1),1);
scatter3(ST2(IR_Count1,1),ST2(IR_Count1,2),ST2(IR_Count1,3),S1,C1(IR_Count1,:),'filled')
xlabel('Tumor cell','FontSize',20,'FontWeight','bold','FontName','Palatino Linotype')
ylabel('Killer T cell','FontSize',20,'FontWeight','bold','FontName','Palatino Linotype')
zlabel('CAF','FontSize',20,'FontWeight','bold','FontName','Palatino Linotype')
xlim([0 1])
ylim([0 1])
zlim([0 1])
set(ax4,'FontName','Palatino Linotype','FontSize',19,'FontWeight','bold');

figure
subplot(1,2,1)
x=1;
bar(x,[median(ST1(IR_Count,1)) median(ST1(IR_Count,2)) median(ST1(IR_Count,3))]);
Q_IR_F_L=quantile(ST1(IR_Count,3),0.25);
Q_IR_F_H=quantile(ST1(IR_Count,3),0.75);
Q_IR_Tum_L=quantile(ST1(IR_Count,1),0.25);
Q_IR_Tum_H=quantile(ST1(IR_Count,1),0.75);
Q_IR_T_L=quantile(ST1(IR_Count,2),0.25);
Q_IR_T_H=quantile(ST1(IR_Count,2),0.75);


subplot(1,2,2)
x=1;
bar(x,[median(ST2(IR_Count1,1)) median(ST2(IR_Count1,2)) median(ST2(IR_Count1,3))]);
Q_IR_F_L_p=quantile(ST2(IR_Count1,3),0.25);
Q_IR_F_H_p=quantile(ST2(IR_Count1,3),0.75);
Q_IR_Tum_L_p=quantile(ST2(IR_Count1,1),0.25);
Q_IR_Tum_H_p=quantile(ST2(IR_Count1,1),0.75);
Q_IR_T_L_p=quantile(ST2(IR_Count1,2),0.25);
Q_IR_T_H_p=quantile(ST2(IR_Count1,2),0.75);