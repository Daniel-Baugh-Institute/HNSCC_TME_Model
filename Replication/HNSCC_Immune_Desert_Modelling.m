%% This code generates the main text plots corresponding to the immune desert TME subtype (Fig. 3 of the manuscript)
%% Initialisation
y_0=[22.0447;6.1926;17.8158;201.5737;5.0551;6.3483;59.8564;7.3397;0;8.5270;4.4444;0.8939;52.8060;15.4442;18.1709;13.8821;4.8854;19.5717;13.3581;10.6108;12.2203;10.1635;16.770;13.7578]; % Initial condition 
alpha=0.005; % Immune-accessibility
P=HNSCC_parameters_immune_desert_modelling(alpha);% Parameters for immune desert TME subtype
Q_ab=P;
Prol_death=[1:2:10]; % Ratio of proliferation to death rate

%% Studying the time profiles for Exhausted and killer T cells pre-ICI for different killer to exhaustion rate (Figure no. 3(a-b))
% Purpose of the study: Colonization of the immune landscape by the exhausted T cells for different exhaustion rate leading to overall immune-desert
% Setting the proliferation rate of the killer-T cells
Q_ab(30)=P(30)*Prol_death(3);
Q_ab(31)=15;
TKPD_TEX=[1:5:35]; % Different conversion rate multiplier from killer to exhausted T cells 
figure
for m=1:1:length(TKPD_TEX)
K_TEX=TKPD_TEX(m)*P(36); % Modifying the killer-exhaustion conversion rate
Q_ab(36)=K_TEX;
subplot(1,2,1)
[t_pre,x_pre]=ode15s(@(t,y)HNSCC_mod(t,y,0,0,0,0,0,0,Q_ab),[0 70000],y_0);
x_pre(x_pre<0)=0;
subplot(1,2,1)
% %yyaxis left
plot(t_pre(1:250),x_pre(1:250,8)/5000,'LineWidth',2);
xlabel('Time','FontSize',18,'FontWeight','bold','FontName','Palatino Linotype')
ylabel('Killer T cells (x_{T_K})','FontSize',18,'FontWeight','bold','FontName','Palatino Linotype')
xlim([0 5]) 
% % yyaxis right 
hold on
subplot(1,2,2)
plot(t_pre,x_pre(:,12)/5000,'LineWidth',2,'LineStyle','-')
xlabel('Time','FontSize',18,'FontWeight','bold','FontName','Palatino Linotype')
ylabel('Exhausted T cells (x_{T_{Exh}})','FontSize',18,'FontWeight','bold','FontName','Palatino Linotype')
hold on
xlim([0 10])
end 
    
%% Study of the T_reg-T_killer trajectory for different CAF-TREG interaction rates       (Figure no. 3(c))
% Purpose of the study: Demonstrating the prospect of a subtype within immune desert with high Regulatory T cell population
Q_C=P;
K_prol=[]; % Proliferation rate of killer T cells
K_death=[]; % Death rate of killer T cells
Prol_death=[1:2:10]; % Ratio of proliferation to death rate
Q_C(30)=P(30)*Prol_death(3); % Setting the proliferation rate of the killer-T cells
Q_C(36)=72;
figure
ax2=axes();
CAF_TREG=[1:2:9];   % CAF-TReg multiplier             
for l=1:1:length(CAF_TREG)
K_CAFTREG=CAF_TREG(l)*P(40); % Modifying the wild type CAF-TREG interaction constant
Q_C(40)=K_CAFTREG;
[t_pre,x_pre]=ode15s(@(t,y)HNSCC_mod(t,y,0,0,0,0,0,0,Q_C),[0 70000],y_0);
x_pre(x_pre<0)=0;
plot(log((1+x_pre(:,11)/5000))/log(2),log2(1+x_pre(:,8)/5000),'LineWidth',2); % The T_Reg population range is too wide to spot the significant transition point with a default representation
xlabel('Regulatory T cells','FontSize',16,'FontWeight','bold','FontName','Palatino Linotype')
ylabel('Killer T cells','FontSize',16,'FontWeight','bold','FontName','Palatino Linotype')
xlim([0 0.0025])
set(ax2,'FontName','Palatino Linotype','FontSize',19,'FontWeight','bold','LineWidth',1);
box off
hold on
end

%% ICI-Response for immune desert: Study of the PD1- killer T-Exhausted T cells for different Exhaustion rates  (Figure No: 3(d))
% Purpose of the study: Demonstrating the efficacy of anti-PD1 in an exhaustion-driven scenario
figure
Prol_death_NPD1=[0.5 5];
Q_d=P;
Q_d(11)=5;
Q_d(16)=1500;
Q_d(30)=Prol_death_NPD1(1)*Q_d(42);
Q_d(31)=Q_d(30);
anti_PD1=1;
[t_post_low_prolif,x_post_low_prolif]=ode23s(@(t,y)HNSCC_mod(t,y,anti_PD1,0,0,0,0,0,Q_d),[0 70000],y_0+200); % For the low-proliferation-driven immune-desert
plot(x_post_low_prolif(:,12)/5000,x_post_low_prolif(:,9)/5000,'LineWidth',2,'LineStyle','--')
xlabel('Exhausted T cells','FontSize',16,'FontWeight','bold','FontName','Palatino Linotype')
ylabel('PD1^- Killer T cells','FontSize',16,'FontWeight','bold','FontName','Palatino Linotype')
hold on
Q_d(30)=Prol_death_NPD1(2)*Q_d(42);
Q_d(31)=Q_d(30);
Q_d(36)=72;
PD1_NPD1=[0.5:5:25];
L=length(PD1_NPD1);
for m=0:1:length(PD1_NPD1)-1
Q_d(35)=PD1_NPD1(L-m);
[t_post,x_post]=ode15s(@(t,y)HNSCC_mod(t,y,anti_PD1,0,0,0,0,0,Q_d),[0 70000],y_0+200);
plot(x_post(:,12)/5000,x_post(:,9)/5000,'LineWidth',2)
xlabel('Exhausted T cells','FontSize',16,'FontWeight','bold','FontName','Palatino Linotype')
ylabel('PD1^- Killer T cells','FontSize',16,'FontWeight','bold','FontName','Palatino Linotype')
hold on
end

%% Study of the PD1- killer T cells and post-ICI tumor cells for different Exhaustion rates (Figure no. 3(e))
% Purpose of the study: Demonstrating the efficacy of anti-PD1 in an exhaustion-driven scenario
figure
Prol_death_NPD1=[0.5 5]; % Proliferation to death balance for Killer T cells
Q_e=P;                   % Storing default parameter values
Q_e(11)=5;               % Adjusting the constant resource rate
Q_e(16)=1500;            % Adjusting the Killer T cell cytotoxicity
Q_e(30)=0.5;             % Low proliferation of killer T cell: Immune-desert due to low proliferation
Q_e(31)=Q_e(30);        
anti_PD1=1;
[t_post_low_prolif,x_post_low_prolif]=ode23s(@(t,y)HNSCC_mod(t,y,anti_PD1,0,0,0,0,0,Q_e),[0 70000],y_0+200);
plot(sum(x_post_low_prolif(:,1:6),2)/30000,x_post_low_prolif(:,9)/5000,'LineWidth',2,'LineStyle','--')
xlabel('Post-ICI tumor cells','FontSize',16,'FontWeight','bold','FontName','Palatino Linotype')  % Ploting Killer T cells vs post-ICI tumor cells in low proliferation-induced immune desert
ylabel('PD1- Killer T cells','FontSize',16,'FontWeight','bold','FontName','Palatino Linotype')
xlim([0 1])
hold on

Q_e(30)=Prol_death_NPD1(2)*Q_e(42);  % Setting up higher proliferation rate
Q_e(31)=Q_e(30);                      
Q_e(36)=72;                          % Ensuring exhaustion-driven immune-desert
PD1_NPD1=[0.5:5:25];                 % Simulating for different anti-PD1 binding rates
L=length(PD1_NPD1);
for m=0:1:length(PD1_NPD1)-1
Q_e(35)=PD1_NPD1(L-m);
[t_post,x_post]=ode15s(@(t,y)HNSCC_mod(t,y,anti_PD1,0,0,0,0,0,Q_e),[0 70000],y_0+200);
plot(sum(x_post(:,1:6),2)/30000,x_post(:,9)/5000,'LineWidth',2)
xlabel('Post-ICI tumor cells','FontSize',16,'FontWeight','bold','FontName','Palatino Linotype') % Ploting Killer T cells vs post-ICI tumor cells in exhaustion-driven immune desert
ylabel('PD1- Killer T cells','FontSize',16,'FontWeight','bold','FontName','Palatino Linotype')
hold on
end

%% In presence of Lactate knockout (Figure no. 3(f))
% Purpose of the study: Lactate improves the overall scenario
figure
C_LAC=10;           % Lactate knockout rate
Prol_death_NPD1=[0.5 5];
Q_f=P;
Q_f(11)=5;
Q_f(16)=1500;
Q_f(30)=Prol_death_NPD1(1)*Q_f(42);
Q_f(31)=Q_f(30);
Q_f(36)=72;
anti_PD1=1;
[t_post_low_prolif,x_post_low_prolif]=ode23s(@(t,y)HNSCC_mod(t,y,anti_PD1,0,0,0,C_LAC,0,Q_f),[0 70000],y_0+200);  % Simulating wihtout lactate
plot(sum(x_post_low_prolif(:,1:6),2)/30000,x_post_low_prolif(:,9)/5000,'LineWidth',2,'LineStyle','--')
xlabel('Post-ICI tumor cells','FontSize',16,'FontWeight','bold','FontName','Palatino Linotype')                   % Ploting Killer T cells vs post-ICI tumor cells in low proliferation-induced immune desert
ylabel('PD1- Killer T cells','FontSize',16,'FontWeight','bold','FontName','Palatino Linotype')
xlim([0 1])
hold on
Q_f(30)=Prol_death_NPD1(2)*Q_f(42);
Q_f(31)=Q_f(30);
Q_f(36)=72;
PD1_NPD1=[0.5:5:25];
L=length(PD1_NPD1);
for m=0:1:length(PD1_NPD1)-1
Q_f(35)=PD1_NPD1(L-m);
[t_post,x_post]=ode15s(@(t,y)HNSCC_mod(t,y,anti_PD1,0,0,0,C_LAC,0,Q_f),[0 70000],y_0+200);
plot(sum(x_post(:,1:6),2)/30000,x_post(:,9)/5000,'LineWidth',2)
xlabel('Post-ICI tumor cells','FontSize',16,'FontWeight','bold','FontName','Palatino Linotype')
ylabel('PD1- Killer T cells','FontSize',16,'FontWeight','bold','FontName','Palatino Linotype')
hold on
end