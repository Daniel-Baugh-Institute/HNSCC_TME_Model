%% This MATLAB code scans through the parameter sets for the HNSCC TME model in a resource-competition setting (Fig S1)
%%% Purpose of the study: The effects of varying resource competition on TME subtypes for varying maximum resource supply rate

%% Initialisation

P_para=2000*lhsdesign(24,10000); % Generating parameter samples
y_0=5*ones(24,1); y_0(9)=0;      % setting up uniform initial conditions
anti_PD1=2;                      % Setting the anti-PD1 level
alpha1=[0.01 0.05 0.1 0.3];      % Varying immune accessibility index
tspan=[0 70000];
CAF_ST=[];                       % Storing pre-ICI steady state CAF population
T_ST=[];                         % Storing pre-ICI steady state Killer T cell population   
TUM_ST=[];                       % Storing pre-ICI steady state total tumor cell population 
n=1;
X_ST=[];                         % Stores Pre-ICI steady states
X_ST_Post=[];                    % Stores Post-ICI steady states
Tum_Prol_D_Bal=[];               % Maximum proliferation-death balance for tumor cells
T_Prol_D_Bal=[];                 % Maximum proliferation-death balance for Killer T cells
CAF_Prol_D_Bal=[];               % Maximum proliferation-death balance CAF
Y_RM=[1 10 10^2 10^3 10^4];      % Maximum resource amounts
for k=1:1:length(alpha1)
    P=HNSCC_parameters(alpha1(k)); 
for g=1:1:length(Y_RM)
    P(8)=0;                       % No spatial competition
    P(11)=Y_RM(g);                % Modulating incoming resource
    P(16)=1500;                   % Default Killer T cell cytotoxicity 
    P(1:6)=P(1:6)/10;             % Low natural proliferation for the tumor cells 
for j=1:1:length(P_para(1,:))
    P(12)=P_para(1,j);            % Tuning Exhausted T cell-induced proliferation of tumor cells
    P(14)=P_para(2,j);            % Tuning CAF-induced proliferation of tumor cells
    P(30)=P_para(3,j);            % Tuning natural proliferation of killer T cells
    P(32)=P_para(4,j);            % Proliferation of Helper T cells
    P(33)=P_para(5,j);            % Proliferation of regulatory T cells
    P(36)=P_para(6,j);            % Killer T cell exhaustion rate
    P(37)=P_para(7,j);            % Helper T cell induced proliferation of killer T cells
    P(38)=P_para(8,j);            % IL-2 induced proliferation rate
    P(39)=P_para(9,j);            % Activation of killer T cells through MHC sensing
    P(40)=P_para(10,j);           % CAF-induced proliferation rate of Tregs
    P(41)=P_para(11,j);           % Tregs induced inhibition of helper T cells
    P(42)=P_para(12,j);           % Death rate of PD1+ killer T cells
    P(43)=P_para(13,j);           % Death rate of helper T cells
    P(44)=P_para(14,j);           % Death rate of PD1+ regulatory T cells
    P(45)=P_para(15,j);           % Death rate of exhausted T cells 
    P(48)=P_para(16,j);           % Proliferation of CAF
    P(49)=P_para(17,j);           % OPN-induced proliferation of CAF 
    P(50)=P_para(18,j);           % Tumor cell-induced growth rate of CAF
    P(51)=P_para(19,j);           % Proximity multiplier between CAF and CAF-protected tumor cells 
    P(53)=P_para(21,j);           % KM for LIF based conversion from wild type to CAF 
    P(54)=P_para(22,j);           % M2-based proliferation of CAF
    P(55)=P_para(23,j);           % LIF-driven conversion rate from FWT to CAFs
    P(58)=P_para(24,j);           % Death rate of CAFs
   [t_pre,x_pre]=ode23s(@(t,y)HNSCC_mod_Res_Comp(t,y,0,0,0,0,0,0,P),tspan,y_0);
   X_ST(j,:)=x_pre(length(x_pre),:);
Tum_Prol_D_Bal(j)=(1+P(12))*(1+P(14))*max(P(1:6))/(P(16)+P(26));
T_Prol_D_Bal(j)=P(30)*(1+P(37))*(1+P(38))/(P(36)+P(42));
CAF_Prol_D_Bal(j)=P(48)*(1+P(49))*(1+P(50))*(1+P(54))/(P(58)+P(56));
end 
% Ensuring the sampled parameters covers all the corner regions in the three dimensional cube (Tum_prol_D_bal, CAF_Prol_D_Bal, T_prol_D_Bal)
for k1=1:1:1000
    P(12)=P_para(1,k1);
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
   X_ST(j+k,:)=x_pre(length(x_pre),:);
   Tum_Prol_D_Bal(j+k1)=P(12)*P(14)*max(P(1:6))/(P(16)+P(26));
   T_Prol_D_Bal(j+k1)=P(30)*P(37)/(P(40)+P(43))*P(38)/(P(36)+P(42));
   CAF_Prol_D_Bal(j+k1)=P(48)*P(49)*P(50)*P(54)/(P(58)+P(56));
end 
Tum_PD=Tum_Prol_D_Bal;
TK_PD=T_Prol_D_Bal;
CAF_PD=CAF_Prol_D_Bal;
Bal_PD=[Tum_PD' TK_PD' CAF_PD'];

C_Tum=sum(X_ST(:,1:6),2)/30000;
C_Tk=sum(X_ST(:,8),2)/5000;
C_CAF=sum(X_ST(:,14),2)/5000;
ST=[C_Tum C_Tk C_CAF];

% Omitting biologically impossible solutions
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


%% Plotting the 3 dimensional scatter plots for tumor cells, CAFS, and killer T cells for different maximum resource capacity
subplot(length(alpha1),length(Y_RM),k*g)
ax3=axes();
C=ST;
%C=[sum(X_ST(1:6))/30000 X_ST(8)/5000 X_ST(14)/5000];
S=20*ones(length(ST),1);
scatter3(log(Bal_PD(:,1)),log(Bal_PD(:,2)),log(Bal_PD(:,3)),S,C,'filled')
xlabel('\zeta_{Tumor cell}','FontSize',20,'FontWeight','bold','FontName','Palatino Linotype')
ylabel('\zeta_{Killer T cell}','FontSize',20,'FontWeight','bold','FontName','Palatino Linotype')
zlabel('\zeta_{CAF}','FontSize',20,'FontWeight','bold','FontName','Palatino Linotype')
set(ax3,'FontName','Palatino Linotype','FontSize',19,'FontWeight','bold');
end
end
