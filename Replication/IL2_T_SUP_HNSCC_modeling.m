%% This code generates the main text figures depicting the effects of IL-2 (Supplementary Figure 6)
%% Initialisation
y_0=[22.0447;6.1926;17.8158;201.5737;5.0551;6.3483;59.8564;1;0;8.5270;4.4444;0.8939;52.8060;15.4442;18.1709;13.8821;0;19.5717;13.3581;10.6108;12.2203;10.1635;16.770;13.7578]; % Initial condition
C_IL2=[0:0.5:2]; 
C_IL21=[2.05:0.02:3];       % IL-2 levels
C_IL2_A=[C_IL2 C_IL21];
alpha=0.005; % Fixing the immune accessibility
P=HNSCC_parameters_IL2(alpha); % Loading the parameters
tspan=[0 70000];
anti_PD1=1;
x_post_CAN=[];
x_post_T=[];
x_post1_CAN=[];
x_post1_T=[];
K_IL2TK=[0.5 1 2 3 4 5];
IL_2U=[0:0.1:1]; % Initial IL-2

for K=1:1:length(K_IL2TK)
P(38)=K_IL2TK(K);
for L=1:1:length(IL_2U)
y_0(17)=IL_2U(L);
[t_post1,x_post1]=ode23s(@(t,y)HNSCC_mod(t,y,anti_PD1,0,0,0,0,0,P),[0 70000],y_0);
subplot(1,length(K_IL2TK),K)
plot(x_post1(:,17)/6000,x_post1(:,9)/5000,'LineWidth',2)                                    % Plotting IL-2 vs. Killer T cells
xlabel('IL-2 level','FontSize',16,'FontWeight','bold','FontName','Palatino Linotype')
ylabel('Killer T cells','FontSize',16,'FontWeight','bold','FontName','Palatino Linotype')
ylim([0 0.0008]); % For the inset
hold on
end
end