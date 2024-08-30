%% This code generates the main text figures depicting the effects of IL-2 (Figure 7)
%% Initialisation
y_0=[22.0447;6.1926;17.8158;201.5737;5.0551;6.3483;59.8564;1;0;8.5270;4.4444;0.8939;52.8060;15.4442;18.1709;13.8821;0;19.5717;13.3581;10.6108;12.2203;10.1635;16.770;13.7578]; % Initial condition
y_1=y_0;
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

%% Killer T cells vs tumor cells for IL-2 intake (Figure No. 7(a))
% Purpose of the study: IL-2 based treatment works in an immune accessible environment
figure
ax2=axes();
for k=1:1:length(C_IL2)
[t_post,x_post]=ode23s(@(t,y)HNSCC_mod(t,y,anti_PD1,C_IL2(k),0,0,0,0,P),[0 70000],y_0);
end

for l=1:1:length(C_IL21)
[t_post1,x_post1]=ode23s(@(t,y)HNSCC_mod(t,y,anti_PD1,C_IL21(l),0,0,0,0,P),[0 70000],y_0);
end
x_post(x_post<0)=0;
x_post1(x_post1<0)=0;
plot(log((x_post1(:,1)+x_post1(:,2)+x_post1(:,3))/(y_0(1)+y_0(2)+y_0(3)))/log(2),log((x_post1(:,9)+x_post1(:,8))/y_0(8))/log(2),'LineWidth',2,'LineStyle','--') % Plotting killer T-Tumor without IL2
xlabel('Total tumor cells','FontSize',16,'FontWeight','bold','FontName','Palatino Linotype')
%ylabel('Killer T cells\\(Without IL2)','FontSize',16,'FontWeight','bold','FontName','Palatino Linotype')
hold on
plot(log((x_post(:,1)+x_post(:,2)+x_post(:,3))/(y_0(1)+y_0(2)+y_0(3)))/log(2),log((x_post(:,9)+x_post(:,8))/y_0(8))/log(2),'LineWidth',2,'LineStyle','-')
%ylabel('Killer T cells\\(With IL2)','FontSize',16,'FontWeight','bold','FontName','Palatino Linotype') % Plotting killer T-Tumor with IL2
ylim([-12 10])
xlim([-12 13])
set(ax2,'FontName','Palatino Linotype','FontSize',19,'FontWeight','bold','LineWidth',1);
box off

%% IL-2 vs Killer T for different initial IL-2 values  (Figure No. 7(b))
% Purpose of the study: A small amount of IL-2 impulse is adequate to displace the TME from an immune deseert to immune hot subtype 
figure
IL_2U=[0:0.1:1]; % Initial IL-2
for L=1:1:length(IL_2U)
    y_0(17)=IL_2U(L);
[t_post1,x_post1]=ode23s(@(t,y)HNSCC_mod(t,y,anti_PD1,0,0,0,0,0,P),[0 70000],y_0);
plot(x_post1(:,17)/6000,x_post1(:,9)/5000,'LineWidth',2)        % Plotting IL-2 vs. Killer T cells
xlabel('IL-2 level','FontSize',16,'FontWeight','bold','FontName','Palatino Linotype')
ylabel('Killer T cells','FontSize',16,'FontWeight','bold','FontName','Palatino Linotype')
ylim([0 0.0008]); % For the inset
hold on
end

%% Killer T cells profiles for different IL-2 intake (Figure No. 7(c))
% Purpose of the study: IL-2-based treatment tap into the positive feedback between IL-2 and Killer T cells
figure
C_IL2=[0:0.5:2];            
C_IL21=[2.05:0.02:3];       % IL-2 levels
C_IL2_A=[C_IL2 C_IL21];
P=HNSCC_parameters_IL2(alpha); % Loading the parameters
for k=1:1:length(C_IL2_A)
[t_post,x_post]=ode23s(@(t,y)HNSCC_mod(t,y,anti_PD1,C_IL2_A(k),0,0,0,0,P),[0 70000],y_1);
x_post(x_post<0)=0;
plot(t_post,(x_post(:,8)+x_post(:,9))/5000,'LineWidth',2) % Plotting killer T cell profiles
xlim([0 15])
xlabel('Time','FontSize',16,'FontWeight','bold','FontName','Palatino Linotype')
hold on
ylabel('Killer T cells','FontSize',16,'FontWeight','bold','FontName','Palatino Linotype')
end
