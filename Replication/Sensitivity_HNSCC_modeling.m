%% This MATLAB code calculates the first-order variance based sensitivity coefficient for the HNSCC TME model
%%% Purpose of the study: Identification of the parameters crucial for tumor cells, killer T cells, and CAF. 

%% Initialisation
num_samples = 5000;                                        % number of parameter samples 
alpha1=0.02;                                               % Operating in moderate to high immune accessibility scenario.
P=HNSCC_parameters(alpha1);                                % Load the nominal parameter values
P=HNSCC_parameters(alpha1);                                % Nominal parameters
P(8)=0;                                                    % No space competition
P(16)=1500;                                                % Toxicity of killer T cells: moderate to high
P(1:6)=P(1:6)/10;                                          % Setting up the nominal natural proliferation rates for the tumor cells
P_para=[];
Para_ID=[12 14 30 32 33 36 37 38 39 40 41 42 43 44 45 48 49 50 51 53 54 55 58]; % Contains the ids of influential parameters
y_0=100*ones(24,1); y_0(9)=0;                              % Uniform initial condition
tspan=[0 70000];                                           % Simulation time
CAF_ST=[];
T_ST=[];
TUM_ST=[];
n=1;
X_ST=[];                                                   % Stores the final population for Tumor cells, killer T cells, and CAF 

%% Generating parameter samples
for k=1:1:length(Para_ID)
    P_para(:,k)=P(Para_ID(k))/3+(3*P(Para_ID(k))-P(Para_ID(k))/3)*lhsdesign(num_samples,1); % Generation of 5000 parameter combinations
end


%% Simulating for different accessibility and parameter sets: Generating outputs for different parameter vectors for unconditional variance
    
    for j=1:1:length(P_para(:,1))
   for k=1:1:length(Para_ID)
       P(Para_ID(k))=P_para(j,k);                                         % Modifying the influential parameters with the parameter samples.  
   end
[t_pre,x_pre]=ode23s(@(t,y)HNSCC_mod(t,y,0,0,0,0,0,0,P),tspan,y_0);
X_ST(j,:)=x_pre(length(x_pre),:);                                         % Storing the final values of population and concentration of all cell states and molecules  
end 

C_Tum=sum(X_ST(:,1:6),2)/30000;                                           % Total normalized tumor cell population
C_Tk=sum(X_ST(:,8),2)/5000;                                               % Total normalized killer T cell population
C_CAF=sum(X_ST(:,14),2)/5000;                                             % Total normalized CAF population    
ST=[C_Tum C_Tk C_CAF];                                                    

%% Conditioning of the final population matrix: eliminating the infeasible solutions

L1=find(ST(:,1)<0);
ST(L1,:)=[];
P_para(L1,:)=[];
L2=find(ST(:,1)>1);
ST(L2,:)=[];
P_para(L2,:)=[];
L3=find(ST(:,2)<0);
ST(L3,:)=[];
P_para(L3,:)=[];
L4=find(ST(:,2)>1);
ST(L4,:)=[];
P_para(L4,:)=[];
L5=find(ST(:,3)<0);
ST(L5,:)=[];
P_para(L5,:)=[];
L6=find(ST(:,3)>1);
ST(L6,:)=[];
P_para(L6,:)=[];

%% Fixing one parameter and varrying the other parameters for conditional variance.

Out_para=[];
X_ST1=[];
for k1=1:1:length(Para_ID)
P_cond=P_para;
C_Tum1=[];                                        
C_Tk1=[];                                         
C_CAF1=[];                                            
ST1=[];   
for j1=1:1:length(P_cond(:,1))
for k2=1:1:length(Para_ID)
P(Para_ID(k2))=P_cond(j1,k);
end
P(Para_ID(k1))=median(P_para(:,k1));
[t_pre,x_pre]=ode23s(@(t,y)HNSCC_mod(t,y,0,0,0,0,0,0,P),tspan,y_0);
X_ST1(j1,:)=x_pre(length(x_pre),:);
end
C_Tum1=sum(X_ST1(:,1:6),2)/30000;                                           % Total normalized tumor cell population
C_Tk1=sum(X_ST1(:,8),2)/5000;                                               % Total normalized killer T cell population
C_CAF1=sum(X_ST1(:,14),2)/5000;                                             % Total normalized CAF population    
ST1=[C_Tum1 C_Tk1 C_CAF1];                                                    
Out_para(:,:,k1)=ST1;
end

%% Calculating sensitivity indices for steady state tumor cell, killer T cells, and CAF population with respect to all the influential parameters 
Var_uncon=var(ST);
eta=[];
ST1=[];
for j1=1:1:length(Para_ID)
ST1=ST;
OUT=[];
OUT=Out_para(:,:,j1);
L1=find(OUT(:,1)<0);
OUT(L1,:)=[];
ST1(L1,:)=[];
L2=find(OUT(:,1)>1);
OUT(L2,:)=[];
ST1(L2,:)=[];
L3=find(OUT(:,2)<0);
OUT(L3,:)=[];
ST1(L3,:)=[];
L4=find(OUT(:,2)>1);
OUT(L4,:)=[];
ST1(L4,:)=[];
L5=find(OUT(:,3)<0);
OUT(L5,:)=[];
ST1(L5,:)=[];
L6=find(OUT(:,3)>1);
OUT(L6,:)=[];
ST1(L6,:)=[];
eta(:,j1)=var(OUT)./var(ST1); % Ratio of the conditional to unconditional variances.
end
