%% This code generates the figure depicting the role of Lactate as biomarker and theraputic target (Figure 10)

%% Initialisation 
alpha=0.005; % Immune accessibility
P=HNSCC_parameters(alpha);
 K_Lac=100;
P(17)=K_Lac;
Q=P;
tspan=[0 70000]
K_TKC=[10 300 600 900 1200 1500 1800 2100 2400 2700 3000 3500 4000 4500 5000]; % Cytotoxic activity for killer T cells
Tum_cell_balance=[]; % Ratio from PDL1-ve to +ve tumor cells
C_Lac1=[0:2:10];
C_Lac2=[15:5:50];
C_Lac=[C_Lac1 C_Lac2];
Lac_ST=[];
%K_Lac=100;

%% Proportion of tumor cells wrt different cytotoxic activity for Lac inhibition rate: Lactate as target (Figure No. 10(a-b))
% 1)Purpose of the study:Lac reprograms the PDL1-/PDL1+ tumor cell ratio towards PDL1- tumor cells

for m=1:1:length(C_LAC)
y_0=[22.0447;6.1926;17.8158;201.5737;5.0551;6.3483;59.8564;7.3397;0;8.5270;4.4444;0.8939;52.8060;15.4442;18.1709;13.8821;4.8854;19.5717;13.3581;10.6108;12.2203;10.1635;16.770;13.7578]; % Initial condition   
for l=1:1:length(K_TKC)
Q(16)=K_TKC(l); % Modifying the cytotoxic activity of killer T cells
[t_pre,x_pre]=ode23s(@(t,y)HNSCC_mod(t,y,0,0,0,0,C_LAC(m),0,Q),tspan,y_0);
y_0=x_pre(length(x_pre(:,1)),:);
Tum_cell_balance(l)=(x_pre(length(x_pre(:,1)),1)+x_pre(length(x_pre(:,1)),3))/x_pre(length(x_pre(:,1)),4); % Calculating the PDL1-/PDL1+ tumor cell proportion
end
plot(K_TKC,Tum_cell_balance,'LineWidth',2)                                     % Plotting tumor cell balance vs. T-cell cytotoxicity levels (Fig. 10(a))
xlabel('T cell cytotoxicity','FontSize',16,'FontWeight','bold','FontName','Palatino Linotype')
ylabel('Tumor cell proportion\\PDL1^{-} vs. PDL1{+}','FontSize',16,'FontWeight','bold','FontName','Palatino Linotype')
hold on
end

C_Lac1=[0:200:900];
C_Lac2=[1500:1000:10^4];
C_LAC=[C_Lac1 C_Lac2];
P=HNSCC_parameters(alpha);
K_TKC=1500;
K_Lac=100;
P(17)=K_Lac;
P(16)=K_TKC;
anti_PD1=1;
figure
ax2=axes();
for m=1:1:length(C_LAC)
y_0=[22.0447;6.1926;17.8158;201.5737;5.0551;6.3483;59.8564;7.3397;0;8.5270;4.4444;0.8939;52.8060;15.4442;18.1709;13.8821;4.8854;19.5717;13.3581;10.6108;12.2203;10.1635;16.770;13.7578]; % Initial condition   
[t_post,x_post]=ode23s(@(t,y)HNSCC_mod(t,y,anti_PD1,0,0,0,C_LAC(m),0,P),tspan,y_0);
plot(t_post,(x_post(:,1)+x_post(:,2)+x_post(:,3)+x_post(:,4)+x_post(:,5)+x_post(:,6))/30000,'LineWidth',2) % Plotting tumor cell balance vs. T-cell cytotoxicity levels
xlabel('Time','FontSize',20.8,'FontWeight','bold','FontName','Palatino Linotype')
ylabel('Tumor cell population','FontSize',20.8,'FontWeight','bold','FontName','Palatino Linotype')
xlim([0 0.03])
ylim([0 1])
hold on
set(ax2,'FontName','Palatino Linotype','FontSize',19,'FontWeight','bold','LineWidth',1);
box off
end


%% Lac vs tumor cells for different accessibility (Figure No. 10(c))
%Purpose of the study: Lactate as a biomarker
figure
y_0=[1089.2;4820.3;3302.5;4106.2;2576.4;0990.7;1340.2;4718.3;0;4673.2;2045.2;684.8;502.5;797.0;2504.0;4912.0;2114.7;1753.9;1601.0;4659.6;4097.4;663.8;3277.8;732.9]; % Initial condition
alpha1=[0.01 0.02 0.04 0.05 0.1 0.3 0.5]; % Fraction of CAFs engadged in building the barrier
Y_LAC=22531; % Maximum LAC Carrying capacity. Can be calculated from the equation describing LIF dynamics
for l=1:1:length(alpha1)
    P=HNSCC_parameters(alpha1(l));
    [t_pre,x_pre]=ode23s(@(t,y)HNSCC_mod(t,y,0,0,0,0,0,0,P),[0 70000],y_0);
    [t_post,x_post]=ode23s(@(t,y)HNSCC_mod(t,y,anti_PD1,0,0,0,0,0,P),[0 70000],x_pre(length(x_pre(:,1)),:));
    I(l)=(1-tanh(alpha1(l)*P(20)*x_pre(length(x_pre(:,1)),14)))*100; %Calculating Immune accessibility index
    Tum_tot=x_post(:,1:6);
     plot(sum(Tum_tot,2)/30000,x_post(:,21)/Y_LAC,'LineWidth',2)                               %Plotting Total tumor cells vs LIF
  xlabel('Tumor cell population','FontSize',16,'FontWeight','bold','FontName','Palatino Linotype')
    ylabel('Lactate concentration','FontSize',16,'FontWeight','bold','FontName','Palatino Linotype')
hold on
end


