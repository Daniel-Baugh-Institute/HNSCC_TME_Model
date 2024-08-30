%% This code generates the main text figures depicting the prospect of IL-8 as biomarker (Figure 9) 
%% Initialisation
y_0=[5.6304;7.6874;8.0805;3.9277;7.5001;0.1805;0.8207;2.5095;0;7.9163;7.5867;8.0238;6.3413;0.4928;9.4814;6.1741;9.7238;7.2642;5.0844;8.1333;7.8140;1.8921;0.9877;5.9869];% initial point

%%  IL-8 vs time (Figure No.- 9(a))
% Purpose of the study: Increase in post-ICI IL-8 can be a potential marker of possible (Non-)response affected by fibro-rich environment
alpha=[0.005 0.08 0.1 0.15 0.2 0.3 0.4 0.5]; % Different accessibility rate
tspan=[0 70000];
I=[];
anti_PD1=2; %ICI
Y_CAN=30000; % Carrying capacity of tumor cells
Y_CAF=5000;  % Carrying capacity of CAF
Y_M2=5000;   % Carrying capacity of M2 macrophage
Y_T=5000;    % Carrying capacity of Killer T cells
IL_8D=5;     % IL-8 degradation rate
IL8_MAX=(2*Y_CAN+2*Y_CAF+15*Y_M2)/IL_8D; % Calculating maximum IL-8 concentration. 2,2,15 are the IL8-secretion rates of tumor cells, CAF, and M2 macrophages respectively
I=[];
for j=1:1:length(alpha)
    P=HNSCC_parameters(alpha(j));
    [t_pre,x_pre]=ode23s(@(t,y)HNSCC_mod(t,y,0,0,0,0,0,0,P),[0 70000],y_0);
    [t_post,x_post]=ode23s(@(t,y)HNSCC_mod(t,y,anti_PD1,0,0,0,0,0,P),[0 70000],x_pre(length(x_pre),:));
    x_post(x_post<0)=0;
    I(j)=(1-tanh(P(20)*alpha(j)*x_pre(length(x_pre),14)))*100; % Calculating immune accessibility
    plot(t_post,x_post(:,20)/IL8_MAX,'LineWidth',2) % plotting IL-8 level profile
     xlabel('Time','FontSize',16,'FontWeight','bold','FontName','Palatino Linotype')
     ylabel('IL-8 level','FontSize',16,'FontWeight','bold','FontName','Palatino Linotype')
     xlim([0 1])
     hold on
end

%% IL-8 for different TME subtypes (Figure No.- 9(b))
% Purpose of the study: IL-8 can be a distinguishing factor between
% different TME subtypes
figure
alpha=[0.005 0.5];
for k=1:1:length(alpha)
P=HNSCC_parameters(alpha(k));
[t_pre,x_pre]=ode23s(@(t,y)HNSCC_mod(t,y,0,0,0,0,0,0,P),[0 70000],y_0);
[t_post,x_post]=ode23s(@(t,y)HNSCC_mod(t,y,anti_PD1,0,0,0,0,0,P),[0 70000],x_pre(length(x_pre),:));
x_post(x_post<0)=0;
plot(x_post(:,9)/Y_T,x_post(:,20)/IL8_MAX,'LineWidth',2) % Plotting IL-8 vs. T cells for immune and fibro rich
xlabel('Killer T cells','FontSize',16,'FontWeight','bold','FontName','Palatino Linotype')
ylabel('IL-8 level','FontSize',16,'FontWeight','bold','FontName','Palatino Linotype')
hold on
end
y_0FD=[2000;0;2500;4000;0;0;59.8564;50;0;40;20;30;150;150;18.1709;13.8821;4.8854;19.5717;13.3581;10.6108;12.2203;10.1635;0;13.7578];
alpha=0.005;
P=HNSCC_parameters_fibro_desert_modelling(alpha);
P(16)=1500;
[t_pre,x_pre]=ode23s(@(t,y)HNSCC_mod(t,y,0,0,0,0,0,0,P),[0 70000],y_0FD);
[t_post,x_post]=ode23s(@(t,y)HNSCC_mod(t,y,anti_PD1,0,0,0,0,0,P),[0 70000],x_pre(length(x_pre),:));
plot(x_post(:,9)/Y_T,x_post(:,20)/IL8_MAX,'LineWidth',2) % Plotting IL8 vs T cells for fibro desert
xlabel('Killer T cells','FontSize',16,'FontWeight','bold','FontName','Palatino Linotype')
ylabel('IL-8 level','FontSize',16,'FontWeight','bold','FontName','Palatino Linotype')
hold on

%% Effect of OPN+LIF knockout on IL-8 levels (Figure No.- 9(c))
figure
y_0=[108.2087;602.6678;350.1559;937.5152;556.5256;33.2473;888.2395;244.3348;0;722.0367;78.8519;697.2312;793.8480;450.4767;541.5014;253.9598;961.2047;294.3169;172.0025;859.1064;766.1546;641.7523;484.2232;132.2313];
C_LIF=[0 500000000];
C_OPN=[0 500000000];
alpha=0.5;
P=HNSCC_parameters_OPN(alpha);
IL_8_Residue=[];
[t_post,x_post]=ode23s(@(t,y)HNSCC_mod(t,y,0,0,C_LIF(1),0,0,C_OPN(1),P),[0 70000],y_0);
y_01=x_post(length(x_post),:);
[t_post_OPN,x_post_OPN]=ode23s(@(t,y)HNSCC_mod(t,y,anti_PD1,0,C_LIF(1),0,0,C_OPN(2),P),[0 70000],y_01);
[t_post_LIF,x_post_LIF]=ode23s(@(t,y)HNSCC_mod(t,y,anti_PD1,0,C_LIF(2),0,0,C_OPN(1),P),[0 70000],y_01);
P(49)=0;
P(55)=0;
[t_post_OPNLIF,x_post_OPNLIF]=ode23s(@(t,y)HNSCC_mod(t,y,anti_PD1,0,C_LIF(2),0,0,C_OPN(1),P),[0 70000],y_01);
x_post_OPNLIF(x_post_OPNLIF<0)=0;
plot(sum(x_post(:,1:6),2)/Y_CAN,x_post(:,20)/IL8_MAX,'LineWidth',2) % Plotting IL-8 vs tumor cells in fibro-rich condition
xlabel('Total tumor cell','FontSize',16,'FontWeight','bold','FontName','Palatino Linotype')
ylabel('IL-8 level','FontSize',16,'FontWeight','bold','FontName','Palatino Linotype')
hold on
plot(sum(x_post_OPN(:,1:6),2)/Y_CAN,x_post_OPN(:,20)/IL8_MAX,'LineWidth',2)% Plotting IL-8 vs tumor cells With LIF knockout
hold on
plot(sum(x_post_LIF(:,1:6),2)/Y_CAN,x_post_LIF(:,20)/IL8_MAX,'LineWidth',2) % Plotting IL-8 vs tumor cells With OPN knockout
hold on
plot(sum(x_post_OPNLIF(:,1:6),2)/Y_CAN,x_post_OPNLIF(:,20)/IL8_MAX,'LineWidth',2) % Plotting IL-8 vs tumor cells With both OPN+LIF knockout
