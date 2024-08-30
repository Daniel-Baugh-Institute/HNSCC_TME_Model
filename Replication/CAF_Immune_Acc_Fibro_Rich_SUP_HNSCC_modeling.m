%% This MATLAB code simulates the model for different levels of immune accessibility (Fig S5)
%%% Purpose of the study: To evaluate how the CAF dynamics modulates the immune accessibility (I_a) over time- sudden/ gradual change in I_a over time  
alpha=0.3;
K_P_CAF1=[0:0.1:1];
K_P_CAF2=[1.1:2:10.1];
K_P_CAF1=[0:0.01:1];
K_P_CAF=[K_P_CAF1 K_P_CAF2];
clc
P=HNSCC_parameters(0.3);
P(11)=20;
y_0=[1089.2;4820.3;3302.5;4106.2;2576.4;0990.7;1340.2;4718.3;0;4673.2;2045.2;684.8;502.5;797.0;2504.0;4912.0;2114.7;1753.9;1601.0;4659.6;4097.4;663.8;3277.8;732.9];
for l=1:1:length(K_P_CAF)
P(49)=K_P_CAF(length(K_P_CAF)-l+1);
P(48)=0.01;
[t_pre,x_pre]=ode23s(@(t,y)HNSCC_mod_Res_Comp(t,y,0,0,0,0,0,0,P),[0 70000],y_0);
I=(1-tanh(alpha*P(20)*x_pre(:,14)));
subplot(1,2,1)
plot(I,x_pre(:,14)/5000,'LineWidth',2)            % Plotting CAF-immune accessibility for different CAF proliferation rate (Fig. S5(a))
xlim([0 1])
ylim([0 1]);
hold on
subplot(1,2,2)
plot(t_pre,I)                                     % Plotting Immune-accessiblity vs time for different CAF proliferation rate (Fig. S5(b))
xlim([0 1])
hold on
end


