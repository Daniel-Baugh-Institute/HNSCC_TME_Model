%% This code generates the main text figures concerning LIF Knockout (Figure S7)

%% Initialisation
y_0=[5.6304;7.6874;8.0805;3.9277;7.5001;0.1805;0.8207;2.5095;0;7.9163;7.5867;8.0238;6.3413;0.4928;9.4814;6.1741;9.7238;7.2642;5.0844;8.1333;7.8140;1.8921;0.9877;5.9869];% Initialisation
tspan=[0 70000];
C_LIF=[0 50 100 150 200 250 300 350 500 650 800 1000 1500 5000 15000000]; % OPN-removal rate
anti_PD1=2; % ICI 
alpha=0.3;% Fraction of CAF blocking T cells

%% Effect of OPN removal (with or without LIF) (Figure No. 8(a-b))
P=HNSCC_parameters_OPN(alpha);
CAN_ACC=[]; % Accessible tumor cells
CAN_IACC=[]; % Inaccessible tumor cells
CAN_ACC_Tot=[]; % Total accessible tumor cells
CAN_IACC_Tot=[];% Total inaccessible tumor cells


for k=1:1:length(C_LIF)
[t_pre,x_pre]=ode23s(@(t,y)HNSCC_mod(t,y,0,0,C_LIF(length(C_LIF)-k+1),0,0,0,P),[0 70000],y_0);
CAN_ACC=[x_pre(:,1) x_pre(:,3:4)];
CAN_IACC=[x_pre(:,2) x_pre(:,5:6)];
CAN_ACC_Tot=sum(CAN_ACC,2);
CAN_IACC_Tot=sum(CAN_IACC,2);
subplot(1,2,1)
plot(CAN_ACC_Tot/30000,CAN_IACC_Tot/30000,'LineWidth',2) % Plotting the propotion of accessible vs inaccessible tumor cells for anti-LIF
xlabel('Accessible tumor cells','FontSize',16,'FontWeight','bold','FontName','Palatino Linotype')
ylabel('Inaccessible tumor cells','FontSize',16,'FontWeight','bold','FontName','Palatino Linotype')
hold on 
end

%% Immune accessibility index wrt OPN levels and decreasing LIF levels (Figure No. 8(c))
% Purpose of the study: OPN removal reprograms the TME towards improved response to ICI
FC_TUM_OPN=[];
FC_TUM_LIF=[];
I_LIF=[];

for k=1:1:length(C_LIF)
[t_pre,x_pre]=ode23s(@(t,y)HNSCC_mod(t,y,0,0,C_LIF(length(C_LIF)-k+1),0,0,0,P),[0 70000],y_0);
%FC_TUM_OPN(j)=sum(x_post(length(x_post(:,1)),1:6))/sum(x_pre(length(x_pre(:,1)),1:6));
I_LIF=1-tanh(P(20)*alpha*x_pre(:,14)); % Calculating the immune accessibility index
subplot(1,2,2)
plot(t_pre, I_LIF,'LineWidth',2) % Plotting I vs time levels for different anti-LIF
 xlabel('Time','FontSize',16,'FontWeight','bold','FontName','Palatino Linotype')
 ylabel('Immune accessibiltiy (I)','FontSize',16,'FontWeight','bold','FontName','Palatino Linotype')
xlim([0 1.5])
 hold on
end
