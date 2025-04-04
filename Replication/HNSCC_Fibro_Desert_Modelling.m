%% This code generates the main text plots corresponding to the fibro desert TME subtype (Figure 4)
%% Initialisation
y_0=[2000;0;2500;4000;0;0;59.8564;50;0;40;20;30;150;150;18.1709;13.8821;4.8854;19.5717;13.3581;10.6108;12.2203;10.1635;0;13.7578]; % initial state for simulation
alpha=0.005; % Immune accessibility index (Not a factor here but for computational issues alpha is kept non-zero)
anti_PD1=2; % ICI concentration
P=HNSCC_parameters_fibro_desert_modelling(alpha); % Parameters for fibro-desert

%% Pre-ICI: Total tumor cells and Killer T cells for different killer-exhaustive conversion to cytotoxic activity ratio. (Figures no. 4(a-b))
% Purpose of the study: The exhausted T cells over all increases the total 
% population of pre-ICI T-cells even in hyper-cytotoxic environment mainly 
% for the increase of PDL1- tumor stem cells.  
Q=P;
K_TexC=[10:6:31]; % Multiplier for the conversion from the killer to exhausted T cells
t_span=[0 70000];
C_total=[]; % Total tumor cells
for k=1:1:length(K_TexC)
K_TXC=P(12);
Q(12)=K_TXC*K_TexC(k); % Modifying the exhausted-tumor interaction rate  
[t_pre,x_pre]=ode15s(@(t,y)HNSCC_mod(t,y,0,0,0,0,0,0,Q),[0 70000],y_0);
C_total=x_pre(:,1:6);
subplot(1,2,1)
plot(t_pre,sum(C_total,2)/30000,'LineWidth',2)                                 % Plotting total tumor cells vs time (Fig. 4(a))   
xlabel('Time','FontSize',16,'FontWeight','bold','FontName','Palatino Linotype')
ylabel('Total tumor cells (x_{Tum})','FontSize',16,'FontWeight','bold','FontName','Palatino Linotype')
xlim([0 1.5])
hold on
subplot(1,2,2)
plot(x_pre(:,8)/5000,(x_pre(:,1)+x_pre(:,3))/20000,'Linewidth',2) % Plotting PDL1- tumor cells vs killer T cells (Fig. 4(b))
xlabel('Killer T cell population (x_{T_K})','FontSize',16,'FontWeight','bold','FontName','Palatino Linotype')
ylabel('Stem and PDL1^{-} tumor cells (x_{St}+x_{PDL1^{-}})','FontSize',16,'FontWeight','bold','FontName','Palatino Linotype')
hold on
end

%% Post ICI: Study of post ICI tumor cell vs T cell population for different constant resources             (Figures no. 4(c-d))
% Purpose of the study: 1)The constant resource rates can lead to recurrence after a succesfull ICI therapy.
%2) ICI reduces the total exhausted T cells count drastically.

Res_Mul=[1:2:10]; % Multiplier for different resource rates
Q(12)=P(12)*K_TexC(2)
figure
for l=1:1:length(Res_Mul)
Q(11)=P(11)*Res_Mul(l); % Modifying Resource rates
[t_post,x_post]=ode15s(@(t,y)HNSCC_mod(t,y,2,0,0,0,0,0,Q),[0 70000],x_pre(length(x_pre),:)); % Post-ICI
subplot(1,3,1)
plot(x_post(:,9)/5000,(x_post(:,1)+x_post(:,4)+x_post(:,3))/30000,'Linewidth',2) %Plotting Post-ICI total tumor cells-Killer T cells trajectory post-ICI (Fig 4(d))
xlabel('Killer T cell population (x_{T_K})','FontSize',16,'FontWeight','bold','FontName','Palatino Linotype')
ylabel('Total tumor cells (x_{Tum})','FontSize',16,'FontWeight','bold','FontName','Palatino Linotype')
ylim([0 0.65])
xlim([0 1])
hold on
subplot(1,3,2)
plot((x_post(:,9)+x_post(:,8))/5000,(x_post(:,12))/5000,'Linewidth',2) % Post-ICI T_exhausted-Killer T trajectory. This plot does not figure in the main manuscript
xlabel('Killer T cell population (x_{T_K})','FontSize',16,'FontWeight','bold','FontName','Palatino Linotype')
ylabel('Exhausted T cells','FontSize',16,'FontWeight','bold','FontName','Palatino Linotype')
hold on
plot((x_pre(:,9)+x_pre(:,8))/5000,(x_pre(:,12))/5000,'Linewidth',2,'LineStyle','--') % Pre-ICI T_exhausted-Killer T trajectory. This plot does not figure in the main manuscript
xlabel('Killer T cell population (x_{T_K})','FontSize',16,'FontWeight','bold','FontName','Palatino Linotype')
ylabel('Exhausted T cells','FontSize',16,'FontWeight','bold','FontName','Palatino Linotype')
subplot(1,3,3)                                                                      
% Ploting the time profiles for killer and exhausted T cells pre- and post- ICI (Fig 4(c))
plot(t_pre(1:150),(x_pre(1:150,9)+x_pre(1:150,8))/5000,'Linewidth',2,'LineStyle','--') %Pre ICI Killer T cells
hold on
plot(t_pre(1:150),x_pre(1:150,12)/5000,'Linewidth',2,'LineStyle','--'); % Pre-ICI Exhausted T cells
t_post1=t_pre(150)+t_post;
hold on
plot(t_post1,(x_post(:,9)+x_post(:,8))/5000,'Linewidth',2) % Post-ICI Killer T cells
hold on
plot(t_post1,x_post(:,12)/5000,'Linewidth',2); % Post-ICI Exhausted T cells
xlabel('Time','FontSize',16,'FontWeight','bold','FontName','Palatino Linotype')
ylabel('T cells population','FontSize',16,'FontWeight','bold','FontName','Palatino Linotype')
xlim([0 1])
end
