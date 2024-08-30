%% This code generates the main text figures concerning the OPN and LIF (Figure 8)

%% Initialisation
y_0=[5.6304;7.6874;8.0805;3.9277;7.5001;0.1805;0.8207;2.5095;0;7.9163;7.5867;8.0238;6.3413;0.4928;9.4814;6.1741;9.7238;7.2642;5.0844;8.1333;7.8140;1.8921;0.9877;5.9869];% Initialisation
tspan=[0 70000];
C_OPN=[100 1000 2000 5000 7000 10000 500000]; % OPN-removal rate
C_LIF=[0 10^8]; % LIF knockout
anti_PD1=2; % ICI 
alpha=0.5;% Fraction of CAF blocking T cells

%% Effect of OPN removal (with or without LIF) (Figure No. 8(a-b))
P=HNSCC_parameters_OPN(alpha);
CAN_ACC=[]; % Accessible tumor cells
CAN_IACC=[]; % Inaccessible tumor cells
CAN_ACC_Tot=[]; % Total accessible tumor cells
CAN_IACC_Tot=[];% Total inaccessible tumor cells
figure;
for k=1:1:length(C_LIF)
for j=1:1:length(C_OPN)
[t_pre,x_pre]=ode23s(@(t,y)HNSCC_mod(t,y,0,0,C_LIF(k),0,0,C_OPN(j),P),[0 70000],y_0);
CAN_ACC=[x_pre(:,1) x_pre(:,3:4)];
CAN_IACC=[x_pre(:,2) x_pre(:,5:6)];
CAN_ACC_Tot=sum(CAN_ACC,2);
CAN_IACC_Tot=sum(CAN_IACC,2);
subplot(1,2,k)
plot(CAN_ACC_Tot/30000,CAN_IACC_Tot/30000,'LineWidth',2) % Plotting the propotion of accessible vs inaccessible tumor cells for OPN levels with/without LIF
xlabel('Accessible tumor cells','FontSize',16,'FontWeight','bold','FontName','Palatino Linotype')
ylabel('Inaccessible tumor cells','FontSize',16,'FontWeight','bold','FontName','Palatino Linotype')
hold on 
end
end
% for the CAF-T barcharts
figure
[t_pre,x_pre]=ode23s(@(t,y)HNSCC_mod(t,y,0,0,C_LIF(1),0,0,C_OPN(1),P),[0 70000],y_0);             % Point p in Fig. 8(a)
CAFT1=[x_pre(length(x_pre),14)/5000 x_pre(length(x_pre),8)/5000];
[t_pre,x_pre]=ode23s(@(t,y)HNSCC_mod(t,y,0,0,C_LIF(1),0,0,C_OPN(length(C_OPN)),P),[0 70000],y_0); % Point q in Fig. 8(a)
CAFT2=[x_pre(length(x_pre),14)/5000 x_pre(length(x_pre),8)/5000];
[t_pre,x_pre]=ode23s(@(t,y)HNSCC_mod(t,y,0,0,C_LIF(2),0,0,C_OPN(length(C_OPN)),P),[0 70000],y_0); % Point r in Fig. 8(b)
CAFT3=[x_pre(length(x_pre),14)/5000 x_pre(length(x_pre),8)/5000];
subplot(1,3,1)
bar(2,CAFT1)
subplot(1,3,2)
bar(2,CAFT2)
subplot(1,3,3)
bar(2,CAFT3)

%% Immune accessibility index wrt OPN levels and decreasing LIF levels (Figure No. 8(c))
% Purpose of the study: OPN removal reprograms the TME towards improved response to ICI
FC_TUM_OPN=[];
FC_TUM_LIF=[];
I_OPN=[];
I_LIF=[];
C_OPN=[500 1000 5000 7000 10200 40000 80000 160000 200000 500000]; % OPN removal rate 
C_LIF=[0 50 100 200 400 500]; % LIF removal rate
Anti_OPN_max=max(C_OPN);
figure
for k=1:1:length(C_LIF)
for j=1:1:length(C_OPN)
[t_pre,x_pre]=ode23s(@(t,y)HNSCC_mod(t,y,0,0,C_LIF(k),0,0,C_OPN(j),P),[0 70000],y_0);
%FC_TUM_OPN(j)=sum(x_post(length(x_post(:,1)),1:6))/sum(x_pre(length(x_pre(:,1)),1:6));
I_OPN(k,j)=1-tanh(P(20)*alpha*x_pre(length(x_pre),14)); % Calculating the immune accessibility index
end
 plot(C_OPN/Anti_OPN_max,I_OPN,'LineWidth',2) % Plotting I vs OPN levels for different LIF
 xlabel('Anti-OPN concentration','FontSize',16,'FontWeight','bold','FontName','Palatino Linotype')
 ylabel('Immune accessibiltiy (I)','FontSize',16,'FontWeight','bold','FontName','Palatino Linotype')
 hold on
end
