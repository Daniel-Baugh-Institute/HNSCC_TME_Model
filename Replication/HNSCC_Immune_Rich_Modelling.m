%%% This code generates the main text plots corresponding to the Immune Rich TME subtype
%% Figure 5 of main manuscript, supplementary figures S2 and S4 
%% Initialisation
y_0=[22.0447;6.1926;17.8158;201.5737;5.0551;6.3483;59.8564;7.3397;0;8.5270;4.4444;0.8939;52.8060;15.4442;18.1709;13.8821;4.8854;19.5717;13.3581;10.6108;12.2203;10.1635;16.770;13.7578]; % Initial condition for Simulation
alpha=0.005; % Immune accessibility
Q=[];
P=HNSCC_parameters(alpha);
P(14)=10;% Tuning CAF-C interactions

Q=P;

%% Pre ICI: CAF vs Killer T cells (Figure No.-5(a), Supplementary figure S2)
%Purpose of the study: 1.) Demonstrating the Proliferation rate of T cells does not necessarily change the steady state population of CAF

T_KProl=[0 0.05 21:10:75]; % Proliferation of killer T cells
L=length(T_KProl);
tspan=[0 70000];
Tk=[];
CAF=[];
Q2=Q;
for k=1:1:length(T_KProl)
Q(30)=T_KProl(k); % Modifying the proliferation rate of T cells
Q2(30)=T_KProl(L-k+1);
[t_pre_1,x_pre_1]=ode15s(@(t,y)HNSCC_mod(t,y,0,0,0,0,0,0,Q2),tspan,y_0);
subplot(1,2,1)
plot(x_pre_1(:,8)/5000,x_pre_1(:,14)/5000,'LineWidth',2)                                            %CAF vs. Killer T cells. Supplementary figure S2
xlabel('PD1^{+} Killer T cells','FontSize',16,'FontWeight','bold','FontName','Palatino Linotype')
ylabel('CAF','FontSize',16,'FontWeight','bold','FontName','Palatino Linotype')
hold on
[t_pre,x_pre]=ode15s(@(t,y)HNSCC_mod(t,y,0,0,0,0,0,0,Q),tspan,y_0);
Tk(k)=x_pre(length(x_pre),8); 
CAF(k)=x_pre(length(x_pre),14);
subplot(1,2,2)
plot((x_pre(:,1)+x_pre(:,3))/30000,x_pre(:,4)/30000,'LineWidth',2)                                   % PDL1+ vs PDL1- tumor cells. (Fig. 5 (a))
xlabel('PDL1^{-} Tumor cells','FontSize',16,'FontWeight','bold','FontName','Palatino Linotype')
ylabel('PDL1^{+} Tumor cells','FontSize',16,'FontWeight','bold','FontName','Palatino Linotype')
hold on
Tk(k)=x_pre(length(x_pre),8);
CAF(k)=x_pre(length(x_pre),14);
end

%% Simulation of helper T cells for different ration between CAF-Tumor and CAF-Treg (Figure No.-5(b))
% Purpose of the study: CAF plays a dual role in regulating the helper T cell population
T_CAF_REG=15;
Q1=P;
Q1(39)=150;
CAF_C_MUL=[0.5 1 1.5 2 2.5 3 4 5 5.5 6 6.5 7 9 11];
TH=[];
Q1(41)=0.001;
figure
Y_0=[654.8;354.5;1227.8;1997.2;860.0;1789.7;29.3;704.7;0;1615.7;1287.4;1536;1395.2;222.9;314;925.6;1869.7;1109.7;521.1;1011.7;1468.4;1724.2;429.1;117];
for k1=1:1:length(CAF_C_MUL)
    Q1(40)=T_CAF_REG;
Q1(14)=T_CAF_REG*CAF_C_MUL(k1);
[t_pre,x_pre]=ode15s(@(t,y)HNSCC_mod(t,y,0,0,0,0,0,0,Q1),tspan,Y_0);
plot(t_pre,x_pre(:,10)/5000,'LineWidth',2)                                                   % Plotting time profile for Helper T cell (Fig. 5(b))
xlim([0 1])
hold on
TH(k1)=x_pre(length(x_pre),10);
end

%% Proportion of tumor cells wrt different cytotoxic activity for different T cell exhaustion rate (Figure No.-5(c))
% 1)Purpose of the study: The exhausted T cells modifies the balance of the PDL1-ve and +ve tumor cells towards PDL1-ve tumor cells

T_KPDexh=[16 32 50 80 128 150]; % Multiplier for conversion from killer T cells to exhausted
K_TKC=[10 300 600 900 1200 1500 1800 2100 2400 2700 3000]; % Cytotoxic activity for killer T cells
Q=P;
Tum_cell_balance=[]; % Ratio from PDL1-ve to +ve tumor cells
figure
for m=1:1:length(T_KPDexh)
    Q(36)=T_KPDexh(m); % Modifying the conversion rate from killer-exhaustion
y_0=[22.0447;6.1926;17.8158;201.5737;5.0551;6.3483;59.8564;7.3397;0;8.5270;4.4444;0.8939;52.8060;15.4442;18.1709;13.8821;4.8854;19.5717;13.3581;10.6108;12.2203;10.1635;16.770;13.7578]; % Initial condition   
for l=1:1:length(K_TKC)
Q(16)=K_TKC(l); % Modifying the cytotoxic activity of killer T cells
[t_pre,x_pre]=ode23s(@(t,y)HNSCC_mod(t,y,0,0,0,0,0,0,Q),tspan,y_0);
y_0=x_pre(length(x_pre(:,1)),:);
Tum_cell_balance(l)=(x_pre(length(x_pre(:,1)),1)+x_pre(length(x_pre(:,1)),3))/x_pre(length(x_pre(:,1)),4);
end
plot(K_TKC,Tum_cell_balance,'LineWidth',2) % Plotting tumor cell balance vs. T-cell cytotoxicity levels
xlabel('PD1^{+} Killer T cells','FontSize',16,'FontWeight','bold','FontName','Palatino Linotype')
ylabel('Tumor cell proportion\\PDL1^{-} vs. PDL1{+}','FontSize',16,'FontWeight','bold','FontName','Palatino Linotype')
hold on
end



%% Post-ICI: Plotting T vs Tumor cells for different CAF-C interaction rate (Figure No.-5(e))
% Purpose of the study: 1) Possibility and elimination of recurrence after successfull ICI
figure
anti_PD1=1;
CAF_C=[0 60 120 180 240 300 500 1000 2000]; % CAF-C interaction rates
P=HNSCC_parameters(0.005);
Q=P;
for n=1:1:length(CAF_C)
    Q(14)=CAF_C(n);% Modifying CAF-C interactions
    [t_pre,x_pre]=ode23s(@(t,y)HNSCC_mod(t,y,0,0,0,0,0,0,Q),tspan,y_0);
    [t_post,x_post]=ode23s(@(t,y)HNSCC_mod(t,y,anti_PD1,0,0,0,0,0,Q),tspan,y_0);
    plot(x_post(:,9)/5000,(x_post(:,1)+x_post(:,3)+x_post(:,4))/30000,'LineWidth',2)% Plotting T vs. total Tumor cells
    xlabel('Killer T cells','FontSize',16,'FontWeight','bold','FontName','Palatino Linotype')
    ylabel('Total tumor cells','FontSize',16,'FontWeight','bold','FontName','Palatino Linotype')
    hold on
end

%% ICI response to immune rich phenotype (Figure No.-5(d))
% Purpose of the study: 1) Is Immune dominated TME responsive to ICI? and
% 2) Is CAF sensitive to ICI in immune dominated region

% Initialisation
num_samples = 10000;                                        % number of parameter samples 
alpha1=0.02;                                               % Operating in moderate to high immune accessibility scenario.
anti_PD1=2;
P=HNSCC_parameters(alpha1);                                % Load the nominal parameter values
P(8)=0;                                                    % No space competition
P(16)=1500;                                                % Toxicity of killer T cells: moderate to high
P(1:6)=P(1:6)/10;                                          % Setting up the nominal natural proliferation rates for the tumor cells
P_para=[];
Para_ID=[12 14 30 32 33 36 37 38 39 40 41 42 43 44 45 48 49 50 51 53 54 55 58]; % Contains the ids of influential parameters
y_0=100*ones(24,1); y_0(9)=0;                              % Uniform initial condition
tspan=[0 70000];                                           % Simulation time
CAF_ST=[];
T_ST=[];
TUM_ST=[];
n=1;
X_ST=[];                                                   % Stores the pre-ICI final population for Tumor cells, killer T cells, and CAF 
X_ST_Post=[];                                              % Stores the pre-ICI final population for Tumor cells, killer T cells, and CAF 
Tum_Prol_D_Bal=[];                                         % Constructing the anti-log(zeta_tumor) vector for low CAF-proliferation region
T_Prol_D_Bal=[];                                           % Constructing the anti-log(zeta_killer T) vector for low CAF-proliferation region
CAF_Prol_D_Bal=[];                                         % Constructing the anti-log(zeta_CAF) vector for low CAF-proliferation region

% Generating parameter samples
for k=1:1:length(Para_ID)
    P_para(:,k)=P(Para_ID(k))/3+(3*P(Para_ID(k))-P(Para_ID(k))/3)*lhsdesign(num_samples,1);      % Generation of 5000 parameter combinations
end    
for j=1:1:length(P_para(:,1))
for k=1:1:length(Para_ID)
     P(Para_ID(k))=P_para(j,k);                                                                  % Modifying the influential parameters with the parameter samples.  
end
[t_pre,x_pre]=ode23s(@(t,y)HNSCC_mod(t,y,0,0,0,0,0,0,P),tspan,y_0);
X_ST(j,:)=x_pre(length(x_pre),:);                                                                % Storing the final values of population and concentration of all cell states and molecules  
[t_post,x_post]=ode23s(@(t,y)HNSCC_mod(t,y,anti_PD1,0,0,0,0,0,P),tspan,x_pre(length(x_pre),:));
X_ST_Post(j,:)=x_post(length(x_post),:);
Tum_Prol_D_Bal(j)=P(12)*P(14)*max(P(1:6))/(P(16)+P(26));                                         % Constructing the anti-log(zeta_tumor) vector for low CAF-proliferation region
T_Prol_D_Bal(j)=P(30)*P(37)/(P(40)+P(43))*P(38)/(P(36)+P(42));                                   % Constructing the anti-log(zeta_killer T) vector for low CAF-proliferation region
CAF_Prol_D_Bal(j)=P(48)*P(49)*P(50)*P(54)/(P(58)+P(56));                                         % Constructing the anti-log(zeta_CAF) vector for low CAF-proliferation region
end 

Tum_PD=Tum_Prol_D_Bal;
TK_PD=T_Prol_D_Bal;
CAF_PD=CAF_Prol_D_Bal;
Bal_PD=[Tum_PD' TK_PD' CAF_PD'];
Bal_PD1=Bal_PD;

C_Tum=sum(X_ST(:,1:6),2)/30000;                                                                  % Total normalized pre_ICI tumor cell population
C_Tk=sum(X_ST(:,8),2)/5000;                                                                      % Total normalized pre_ICI killer T cell population
C_CAF=sum(X_ST(:,14),2)/5000;                                                                    % Total normalized pre_ICI CAF population    
ST=[C_Tum C_Tk C_CAF];                                                    

C_Tum_Post=sum(X_ST_Post(:,1:6),2)/30000;                                                        % Total normalized post_ICI tumor cell population
C_Tk_Post=sum(X_ST_Post(:,8),2)/5000;                                                            % Total normalized post_ICI killer T cell population
C_CAF_Post=sum(X_ST_Post(:,14),2)/5000;                                                          % Total normalized post_ICI CAF population    
ST_Post=[C_Tum_Post C_Tk_Post C_CAF_Post];                                                    

% Conditioning of the final population matrix: eliminating the infeasible solutions

L1=find(ST(:,1)<0);
ST(L1,:)=[];
Bal_PD1(L1,:)=[];
L2=find(ST(:,1)>1);
ST(L2,:)=[];
Bal_PD1(L2,:)=[];
L3=find(ST(:,2)<0);
ST(L3,:)=[];
Bal_PD1(L3,:)=[];
L4=find(ST(:,2)>1);
ST(L4,:)=[];
Bal_PD1(L4,:)=[];
L5=find(ST(:,3)<0);
ST(L5,:)=[];
Bal_PD1(L5,:)=[];
L6=find(ST(:,3)>1);
ST(L6,:)=[];
Bal_PD1(L6,:)=[];

L1_Post=find(ST_Post(:,1)<0);
ST_Post(L1_Post,:)=[];
Bal_PD(L1_Post,:)=[];
L2_Post=find(ST_Post(:,1)>1);
ST_Post(L2_Post,:)=[];
Bal_PD(L2_Post,:)=[];
L3_Post=find(ST_Post(:,2)<0);
ST_Post(L3_Post,:)=[];
Bal_PD(L3_Post,:)=[];
L4_Post=find(ST_Post(:,2)>1);
ST_Post(L4_Post,:)=[];
Bal_PD(L4_Post,:)=[];
L5_Post=find(ST_Post(:,3)<0);
ST_Post(L5_Post,:)=[];
Bal_PD(L5_Post,:)=[];
L6_Post=find(ST_Post(:,3)>1);
ST_Post(L6_Post,:)=[];
Bal_PD(L6_Post,:)=[];

% Plotting the zeta values for tumor cells, Killer T cells, and CAF population Pre-ICI: 
% Each point is color-coded by the respective tumor cells, killer T cells, and CAF population
figure
ax3=axes();
C=ST;
S=20*ones(length(ST),1);
scatter3(log(Bal_PD1(:,1)),log(Bal_PD1(:,2)),log(Bal_PD1(:,3)),S,C,'filled')
xlabel('\zeta_{Tumor cell}','FontSize',20,'FontWeight','bold','FontName','Palatino Linotype')
ylabel('\zeta_{Killer T cell}','FontSize',20,'FontWeight','bold','FontName','Palatino Linotype')
zlabel('\zeta_{CAF}','FontSize',20,'FontWeight','bold','FontName','Palatino Linotype')
set(ax3,'FontName','Palatino Linotype','FontSize',19,'FontWeight','bold');

% Plotting the zeta values for tumor cells, Killer T cells, and CAF population Post-ICI: 
% Each point is color-coded by the respective tumor cells, killer T cells, and CAF population
figure
ax3=axes();
C=ST_Post;
S=20*ones(length(ST_Post),1);
scatter3(log(Bal_PD(:,1)),log(Bal_PD(:,2)),log(Bal_PD(:,3)),S,C,'filled')
xlabel('\zeta_{Tumor cell}','FontSize',20,'FontWeight','bold','FontName','Palatino Linotype')
ylabel('\zeta_{Killer T cell}','FontSize',20,'FontWeight','bold','FontName','Palatino Linotype')
zlabel('\zeta_{CAF}','FontSize',20,'FontWeight','bold','FontName','Palatino Linotype')
set(ax3,'FontName','Palatino Linotype','FontSize',19,'FontWeight','bold');

%% Plotting CAF-FWT population post ICI (Figure No.-Supplementary figure 4)
% Purpose of study: Immune-Rich TME reverse the pre-ICI CAF-FWT balance: A characteristics of Immune-Rich subtype 
Q=P;
anti_PD1=1;
tspan=[0 70000];
y_0=[22.0447;6.1926;17.8158;201.5737;5.0551;6.3483;59.8564;7.3397;0;8.5270;4.4444;0.8939;52.8060;15.4442;18.1709;13.8821;4.8854;19.5717;13.3581;10.6108;12.2203;10.1635;16.770;13.7578]; % Initial condition for Simulation
[t_pre,x_pre]=ode23s(@(t,y)HNSCC_mod(t,y,0,0,0,0,0,0,Q),tspan,y_0);
[t_post,x_post]=ode23s(@(t,y)HNSCC_mod(t,y,anti_PD1,0,0,0,0,0,Q),tspan,x_pre(length(x_pre),:));
%yyaxis left
figure
ax2=axes();
L1=max(find(t_pre<=0.8));
t_post_1=t_post+t_pre(L1);
plot(t_pre(1:L1),x_pre(1:L1,14)/5000,'LineWidth',2,'LineStyle','--') % Plotting pre-ICI CAF population wrt time
hold on
plot(t_post_1,x_post(:,14)/5000,'LineWidth',2,'LineStyle','-') % Plotting post-ICI CAF population wrt time
xlim([0 1.5])
xlabel('Time','FontSize',16,'FontWeight','bold','FontName','Palatino Linotype')
ylabel('CAF population','FontSize',16,'FontWeight','bold','FontName','Palatino Linotype')
set(ax2,'FontName','Palatino Linotype','FontSize',19,'FontWeight','bold','LineWidth',1);
%yyaxis right
figure
ax3=axes();
plot(t_pre(1:L1),x_pre(1:L1,13)/5000,'LineWidth',2,'LineStyle','--') % Plotting pre-ICI Wild type fibroblasts population wrt time
hold on
plot(t_post_1,x_post(:,13)/5000,'LineWidth',2,'LineStyle','-') % Plotting post-ICI Wild type fibroblasts population wrt time
xlim([0 1.5])
xlabel('Time','FontSize',16,'FontWeight','bold','FontName','Palatino Linotype')
ylabel('Tumor suppressive fibroblasts','FontSize',16,'FontWeight','bold','FontName','Palatino Linotype')
set(ax3,'FontName','Palatino Linotype','FontSize',19,'FontWeight','bold','LineWidth',1);