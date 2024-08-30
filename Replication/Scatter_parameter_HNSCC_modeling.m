%% This MATLAB code scans through the parameter sets for the HNSCC TME model: Fig 2 of the manuscript 
%%% Purpose of the study: Identify possible pre-ICI cohorts that can
%%% explain the distinct clinically observed TME subtypes
%%% Assumption: Constant resource supply=> No resource competition

%% Initialisation
P_para=2000*lhsdesign(24,10000);                           % Generate a sample of 10000 different parameter vectors of dimension 24 from 0-2000.
y_0=5*ones(24,1); y_0(9)=0;                                % Uniform initial condition
%anti_PD1=2;                                                % Uniform ICI dosage
alpha1=0.02;                                               % Fixing the immune accessibility 
tspan=[0 70000];                                           % Simulation time
CAF_ST=[];
T_ST=[];
TUM_ST=[];
n=1;
X_ST=[];                                                   % Stores the final population for Tumor cells, killer T cells, and CAF 
X_ST_Post=[];
Tum_Prol_D_Bal=[];                                         % Stores zeta_Tumor cells
T_Prol_D_Bal=[];                                           % Stores zeta_Killer_T
CAF_Prol_D_Bal=[];                                         % Stores zeta_CAF

%% Simulating for different accessibility and parameter sets
for k=1:1:length(alpha1)
    P=HNSCC_parameters(alpha1(k));                         % Nominal parameters
    P(8)=0;                                                % No space competition
    P(16)=1500;                                            % Toxicity of killer T cells: moderate to high
    P(1:6)=P(1:6)/10;                                      % Setting up the nominal natural proliferation rates for the tumor cells
    for j=1:1:length(P_para(1,:))
    P(12)=P_para(1,j);                                     % Rate of exhausted T cell-driven growth of tumor cells 
    P(14)=P_para(2,j);                                     % CAF-driven growth of tumor cells

    P(30)=P_para(3,j);                                     % Proliferation rate of Killer T cells 
    P(32)=P_para(4,j);                                     % Proliferation rate of helper T cells
    P(33)=P_para(5,j);                                     % Proliferation rate of regulator T cells   
    P(36)=P_para(6,j);                                     % Conversion rate from PD1+ T cels to exhausted T cells
    P(37)=P_para(7,j);                                     % Helper-driven proliferation rate of Killer T cells
    P(38)=P_para(8,j);                                     % IL-2 induced proliferation rate of killer T cells
    P(39)=P_para(9,j);                                     % Tumor cell-induced activation of helper T cells via antigen sensing
    P(40)=P_para(10,j);                                    % CAF induced proliferation of regulatory T cells 
    P(41)=P_para(11,j);                                    % TREG induced inhibition of helper T cells  
    P(42)=P_para(12,j);                                    % Death rate of killer PD1+ T cells
    P(43)=P_para(13,j);                                    % Death rate of helper T cells
    P(44)=P_para(14,j);                                    % Death rate of regulatory T cells
    P(45)=P_para(15,j);                                    % Death rate of exhausted T cells

    P(48)=P_para(16,j);                                    % Proliferation rate of CAF cells
    P(49)=P_para(17,j);                                    % OPN-induced Proliferation rate of CAF cells
    P(50)=P_para(18,j);                                    % Tumor cells-induced proliferation rate of CAF
    P(51)=P_para(19,j);                                    % Proximity factor for CAF and CAF-protected tumor cells
    P(53)=P_para(21,j);                                    % Hill like dissociation constant for LIF-driven conversion from CAF to FWT
    P(54)=P_para(22,j);                                    % M2-driven proliferation rate for CAFs 
    P(55)=P_para(23,j);                                    % LIF-driven conversion rate from FWT to CAFs
    P(58)=P_para(24,j);                                    % Death rate for CAF  
   
[t_pre,x_pre]=ode23s(@(t,y)HNSCC_mod(t,y,0,0,0,0,0,0,P),tspan,y_0);

X_ST(j,:)=x_pre(length(x_pre),:);                                         % Storing the final values of population and concentration of all cell states and molecules  
Tum_Prol_D_Bal(j)=(1+P(12))*(1+P(14))*max(P(1:6))/(P(16)+P(26));          % Constructing the anti-log(zeta_tumor) vector
T_Prol_D_Bal(j)=P(30)*(1+P(37))*(1+P(38))/(P(36)+P(42));                  % Constructing the anti-log(zeta_killer T) vector
CAF_Prol_D_Bal(j)=P(48)*(1+P(49))*(1+P(50))*(1+P(54))/(P(58)+P(56));      % Constructing the anti-log(zeta_CAF) vector
end 
end

%%% Simulating for different parameter sets pertaining to fibro desert and desert
%%% Purpose of the study: The initial parameter sampling does niot capture enough points for low over all CAF proliferation rate
% Using the same parameter values for all parameters other than the ones explicitly connected to the dynamics of CAF

for k=1:1:1000
    P(12)=P_para(1,k);
    P(30)=P_para(2,k);
    P(32)=P_para(3,k);
    P(33)=P_para(4,k);
    P(36)=P_para(5,k);
    P(37)=P_para(6,k);
    P(38)=P_para(7,k);
    P(39)=P_para(8,k);
    P(40)=P_para(9,k);
    P(41)=P_para(10,k);
    P(42)=P_para(11,k);
    P(43)=P_para(12,k);
    P(44)=P_para(13,k);
    P(45)=P_para(14,k);
    P(58)=P_para(15,k); 
    P(47)=P(58)*lhsdesign(1,1);                                           % Ensuring fibro-desert region
    P(48)=10*lhsdesign(1,1);
    P(49)=10*lhsdesign(1,1);
    P(50)=10*lhsdesign(1,1);
    P(51)=10*lhsdesign(1,1);
    P(52)=P_para(19,k);
    P(53)=P_para(20,k);
    P(54)=10*lhsdesign(1,1);
    P(55)=10*lhsdesign(1,1);
   [t_pre,x_pre]=ode23s(@(t,y)HNSCC_mod(t,y,0,0,0,0,0,0,P),tspan,y_0);
   X_ST(j+k,:)=x_pre(length(x_pre),:);                                         
   Tum_Prol_D_Bal(j+k)=P(12)*P(14)*max(P(1:6))/(P(16)+P(26));             % Constructing the anti-log(zeta_tumor) vector for low CAF-proliferation region
   T_Prol_D_Bal(j+k)=P(30)*P(37)/(P(40)+P(43))*P(38)/(P(36)+P(42));       % Constructing the anti-log(zeta_killer T) vector for low CAF-proliferation region
   CAF_Prol_D_Bal(j+k)=P(48)*P(49)*P(50)*P(54)/(P(58)+P(56));             % Constructing the anti-log(zeta_CAF) vector for low CAF-proliferation region
end 

Tum_PD=Tum_Prol_D_Bal;
TK_PD=T_Prol_D_Bal;
CAF_PD=CAF_Prol_D_Bal;
Bal_PD=[Tum_PD' TK_PD' CAF_PD'];

C_Tum=sum(X_ST(:,1:6),2)/30000;                                           % Total normalized tumor cell population
C_Tk=sum(X_ST(:,8),2)/5000;                                               % Total normalized killer T cell population
C_CAF=sum(X_ST(:,14),2)/5000;                                             % Total normalized CAF population    
ST=[C_Tum C_Tk C_CAF];                                                    
L1=find(ST(:,1)<0);

%% Conditioning of the final population matrix: eliminating the infeasible solutions
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


%% Plotting the zeta values for tumor cells, Killer T cells, and CAF population: 
% Each point is color-coded by the respective tumor cells, killer T cells, and CAF population
figure
ax3=axes();
C=ST;
S=20*ones(length(ST),1);
scatter3(log(Bal_PD(:,1)),log(Bal_PD(:,2)),log(Bal_PD(:,3)),S,C,'filled')
xlabel('\zeta_{Tumor cell}','FontSize',20,'FontWeight','bold','FontName','Palatino Linotype')
ylabel('\zeta_{Killer T cell}','FontSize',20,'FontWeight','bold','FontName','Palatino Linotype')
zlabel('\zeta_{CAF}','FontSize',20,'FontWeight','bold','FontName','Palatino Linotype')
set(ax3,'FontName','Palatino Linotype','FontSize',19,'FontWeight','bold');


%% Plotting the barcharts showing respective cohorts.

% Fibro-dominated/ rich
FR_Count=find(ST(:,1)>0.05 & ST(:,2)>0.05 & ST(:,3)>0.2);
x=1;
bar(x,[mean(ST(FR_Count,1)) mean(ST(FR_Count,2)) mean(ST(FR_Count,3))]);
Q_FR_F_L=quantile(ST(FR_Count,3),0.25);
Q_FR_F_H=quantile(ST(FR_Count,3),0.75);
Q_FR_Tum_L=quantile(ST(FR_Count,1),0.25);
Q_FR_Tum_H=quantile(ST(FR_Count,1),0.75);
Q_FR_T_L=quantile(ST(FR_Count,2),0.25);
Q_FR_T_H=quantile(ST(FR_Count,2),0.75);

% Immune-dominated/ rich
IR_Count=find(ST(:,2)>0.05 & ST(:,1)<0.6 & ST(:,3)>0.2);
x=1;
bar(x,[mean(ST(IR_Count,1)) mean(ST(IR_Count,2)) mean(ST(IR_Count,3))]);
Q_IR_F_L=quantile(ST(IR_Count,3),0.25);
Q_IR_F_H=quantile(ST(IR_Count,3),0.75);
Q_IR_Tum_L=quantile(ST(IR_Count,1),0.25);
Q_IR_Tum_H=quantile(ST(IR_Count,1),0.75);
Q_IR_T_L=quantile(ST(IR_Count,2),0.25);
Q_IR_T_H=quantile(ST(IR_Count,2),0.75);

% Fibro-desert
FD_Count=find(ST(:,3)<0.005 & ST(:,2)>0.1);
x=1;
bar(x,[mean(ST(FD_Count,1)) mean(ST(FD_Count,2)) mean(ST(FD_Count,3))]);
Q_FD_F_L=quantile(ST(FD_Count,3),0.25);
Q_FD_F_H=quantile(ST(FD_Count,3),0.75);
Q_FD_Tum_L=quantile(ST(FD_Count,1),0.25);
Q_FD_Tum_H=quantile(ST(FD_Count,1),0.75);
Q_FD_T_L=quantile(ST(FD_Count,2),0.25);
Q_FD_T_H=quantile(ST(FD_Count,2),0.75);

% Desert
D_Count=find(ST(:,2)<0.1&ST(:,3)<0.1);
x=1;
bar(x,[mean(ST(D_Count,1)) mean(ST(D_Count,2)) mean(ST(D_Count,3))]);
Q_D_F_L=quantile(ST(D_Count,3),0.25);
Q_D_F_H=quantile(ST(D_Count,3),0.75);
Q_D_Tum_L=quantile(ST(D_Count,1),0.25);
Q_D_Tum_H=quantile(ST(D_Count,1),0.75);
Q_D_T_L=quantile(ST(D_Count,2),0.25);
Q_D_T_H=quantile(ST(D_Count,2),0.75);
