%% This code generates the main text figures corresponding to the Fibro-Rich TME subtype (Figure 6)
%% Initialisation
y_0=[1089.2;4820.3;3302.5;4106.2;2576.4;0990.7;1340.2;4718.3;0;4673.2;2045.2;684.8;502.5;797.0;2504.0;4912.0;2114.7;1753.9;1601.0;4659.6;4097.4;663.8;3277.8;732.9];% Initial condition
tspan=[0 70000];

%% Tumor cells fold change for different immune accessibility indices (Figure No. 6(b))
% Purpose of the study: To illusttrate the effect of the immune accessibility on the overall prognostics for different Cytotoxicity
alpha=[0.01 0.02 0.05 0.08 0.09 0.095 0.1 0.15 0.2 0.25 0.3 0.5]; % fraction of CAF engaged in bilding the barrier
T_KC=[50 100 300 600 900 1200 1500]; % Cytotoxicity of the killer T cells
anti_PD1=1; %ICI
Q=[];
FC_TUM=[]; % Tumor fold change pre vs. post ICI
I=[]; % Immune accessibility index
Lac=[]; % Lactate
y_01=y_0;
P=HNSCC_parameters(0.05);
rho_Cyto_Tex=T_KC/P(12); % Ratio of cytotoxic to oncogenic (Exhausted T-> Tumor cells) activity of T cells 
for k=1:1:length(T_KC)
    for j=1:1:length(alpha)
    P=HNSCC_parameters(alpha(j));
    Q=P;
    Q(16)=T_KC(k); % Modifying the cytotoxic activity
    [t_pre,x_pre]=ode23s(@(t,y)HNSCC_mod(t,y,0,0,0,0,0,0,Q),[0 70000],y_0);
    [t_post,x_post]=ode23s(@(t,y)HNSCC_mod(t,y,anti_PD1,0,0,0,0,0,Q),[0 70000],x_pre(length(x_pre(:,1)),:));
    FC_TUM(j)=sum(x_post(length(x_post),1:6))/sum(x_pre(length(x_pre),1:6)); % Calculating the fold change
    I(j)=(1-tanh(alpha(j)*Q(20)*x_pre(length(x_pre),14)))*100; % Calculating immune accessibility index
    y_0=x_pre(length(x_pre),:);
    y_01=x_post(length(x_post),:);
    end
    plot(I,FC_TUM,'LineWidth',2) % Plotting tumor fold change for different immune accessibility
    hold on
end

%% LIF vs tumor cells for different accessibility (Figure No. 6(c))   
%Purpose of the study: To demonstrate the reduction/stagnation of LIF plays the governing role in the post-ICI CAF population
figure
y_0=[1089.2;4820.3;3302.5;4106.2;2576.4;0990.7;1340.2;4718.3;0;4673.2;2045.2;684.8;502.5;797.0;2504.0;4912.0;2114.7;1753.9;1601.0;4659.6;4097.4;663.8;3277.8;732.9]; % Initial condition
alpha1=[0.01 0.02 0.04 0.05 0.1 0.3 0.5]; % Fraction of CAFs engadged in building the barrier
Y_LIF=22531; % Maximum LIF Carrying capacity. Can be calculated from the equation describing LIF dynamics
for l=1:1:length(alpha1)
    P=HNSCC_parameters(alpha1(l)); % Simulating for different immune accessibility indices
    [t_pre,x_pre]=ode23s(@(t,y)HNSCC_mod(t,y,0,0,0,0,0,0,P),[0 70000],y_0);
    [t_post,x_post]=ode23s(@(t,y)HNSCC_mod(t,y,anti_PD1,0,0,0,0,0,P),[0 70000],x_pre(length(x_pre(:,1)),:));
    I(l)=(1-tanh(alpha1(l)*P(20)*x_pre(length(x_pre(:,1)),14)))*100; %Calculating Immune accessibility index
    Tum_tot=x_post(:,1:6);  % Total post-ICI tumor cells
     plot(x_post(:,18)/Y_LIF,sum(Tum_tot,2)/30000,'LineWidth',2)               %Plotting post-ICI Total tumor cells vs LIF
  xlabel('Total tumor cells','FontSize',16,'FontWeight','bold','FontName','Palatino Linotype')
    ylabel('LIF','FontSize',16,'FontWeight','bold','FontName','Palatino Linotype')
hold on
end

%% CAF profiles for different accesssibility (Figure No. 6(d))            
%Purpose of the study: Changing CAF profile vis-a-vis pre and post ICI is a key differentiator between immune and fibro rich 
figure
y_0=[1089.2;4820.3;3302.5;4106.2;2576.4;0990.7;1340.2;4718.3;0;4673.2;2045.2;684.8;502.5;797.0;2504.0;4912.0;2114.7;1753.9;1601.0;4659.6;4097.4;663.8;3277.8;732.9]; % Initial condition
alpha1=[0.01 0.02 0.04 0.05 0.1 0.3 0.5]; % Fraction of CAFs engadged in building the barrier
for l=1:1:length(alpha1)
    P=HNSCC_parameters(alpha1(l));
    [t_pre,x_pre]=ode23s(@(t,y)HNSCC_mod(t,y,0,0,0,0,0,0,P),[0 70000],y_0);
    [t_post,x_post]=ode23s(@(t,y)HNSCC_mod(t,y,anti_PD1,0,0,0,0,0,P),[0 70000],x_pre(length(x_pre(:,1)),:));
    y_0=x_pre(length(x_pre),:);
    I(l)=(1-tanh(alpha1(l)*P(20)*x_pre(length(x_pre(:,1)),14)))*100; % Calculating immune accessibility index
plot(t_post,x_post(:,14)/5000,'LineWidth',2)                            % Plotting post-ICI CAF vs time
xlabel('Time','FontSize',16,'FontWeight','bold','FontName','Palatino Linotype')
ylabel('Post-ICI CAF','FontSize',16,'FontWeight','bold','FontName','Palatino Linotype')
xlim([0 5])
hold on
end


