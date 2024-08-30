%% This code generates the Supplementary plots corresponding to the Immune Rich TME subtype (S2 and S3) in resource-competition setting
%% Initialisation
y_0=[22.0447;6.1926;17.8158;201.5737;5.0551;6.3483;59.8564;7.3397;0;8.5270;4.4444;0.8939;52.8060;15.4442;18.1709;13.8821;4.8854;19.5717;13.3581;10.6108;12.2203;10.1635;16.770;13.7578]; % Initial condition for Simulation
alpha=0.005; % Immune accessibility
Q=[];
P=HNSCC_parameters(alpha);
P(14)=10;% Tuning CAF-C interactions
P(11)=10^2;
Q=P;

%% Pre ICI: CAF vs Killer T cells 
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
plot(x_pre_1(:,8)/5000,x_pre_1(:,14)/5000,'LineWidth',2)                                  %CAF vs. Killer T cells (Fig. S2)
xlabel('PD1^{+} Killer T cells','FontSize',16,'FontWeight','bold','FontName','Palatino Linotype')
ylabel('CAF','FontSize',16,'FontWeight','bold','FontName','Palatino Linotype')
hold on
[t_pre,x_pre]=ode15s(@(t,y)HNSCC_mod(t,y,0,0,0,0,0,0,Q),tspan,y_0);
Tk(k)=x_pre(length(x_pre),8); 
CAF(k)=x_pre(length(x_pre),14);
Tk(k)=x_pre(length(x_pre),8);
CAF(k)=x_pre(length(x_pre),14);
end

%% Simulation of helper T cells for different ration between CAF-Tumor and CAF-Treg for different resource rates (Figure S3)
% Purpose of the study: CAF plays a dual role in regulating the helper T cell population
T_CAF_REG=15;
Res=[1:2:10];
CAF_C_MUL=[0.5 1 1.5 2 2.5 3 4 5 5.5 6 6.5 7 9 11];
TH=[];
Y_0=[654.8;354.5;1227.8;1997.2;860.0;1789.7;29.3;704.7;0;1615.7;1287.4;1536;1395.2;222.9;314;925.6;1869.7;1109.7;521.1;1011.7;1468.4;1724.2;429.1;117];
for k=1:1:length(Res)
P(11)=Res(k);
Q1=P;
Q1(39)=150;
Q1(41)=0.001;

for k1=1:1:length(CAF_C_MUL)
 Q1(40)=T_CAF_REG;
Q1(14)=T_CAF_REG*CAF_C_MUL(k1);
[t_pre,x_pre]=ode15s(@(t,y)HNSCC_mod(t,y,0,0,0,0,0,0,Q1),tspan,Y_0);
subplot(1,length(Res),k)
plot(t_pre,x_pre(:,10)/5000,'LineWidth',2)   % Plotting Helper T cells vs time          
xlim([0 1])
hold on
TH(k1)=x_pre(length(x_pre),10);
end
end