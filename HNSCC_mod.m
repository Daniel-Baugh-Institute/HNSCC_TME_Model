function dydt=HNSCC_mod(t,y,u,C_IL2,C_LIF,C_IL8,C_Lac,C_OPN,P)

%% Load the parameters

%% Resource to cancer cells

K_RCST=P(1);               %  per unit of resource per unit of time. Resource to T-exposed tumor stem cell 

K_RCSNT=P(2);              %  per unit of resource per unit of time. Resource to Non-T-exposed tumor stem cell

K_RCNPDL1=P(3);            %  per unit of resource per unit of time. Resource to T-exposed PDL1- tumor cell 

K_RCRNPDL1=P(4);           %  per unit of resource per unit of time. Resource to Non-T-exposed PDL1- tumor cell

K_RCPDL1=P(5);             %  per unit of resource per unit of time. Resource to T-exposed PDL1+ tumor cell

K_RCRPDL1=P(6);            %  per unit of resource per unit of time. Resource to Non-T-exposed PDL1+ tumor cell

K_RIn=P(7);                %  Resource per unit time. Default resource intake rate  

alpha_Comp=P(8); 
%% Cancer cell proliferation, interactions, and death

Y_CST=P(9);            % Carrying capacity for PD1- T cells

Y_CPDL1M=P(10);             % Carrying capacity for PD1+ T cells

Y_RM=P(11);                % Resources. Maximum resource withholding capacity.

K_TXC=P(12);               % Exhausted T cells-driven growth of tumor cells.

K_FWTC=P(13);              % Fraction of wild type fibroblasts in proximity with the tumor cells.

K_CAFC= P(14);             % CAF-driven growth of the tumor cells.

K_CAFCR=P(15);             % CAF-driven growth of the tumor non-T exposed cells.
            
K_TKC=P(16);               % T cells-driven apoptosis of tumor cells.

K_Lac=P(17);               % Inhibition rate of lactate

K_TKCIFNG=P(18);           % IFNG induced T cells-driven apoptosis of tumor cells.

alpha=P(19);               % CAF barrier

K_CAFB=P(20);              % Barrier formation rate

K_CSTCNPDL1=P(21);         % Stem to PDL1- 

K_CSTCPDL1=P(22);          % Stem to PDL1+ 

K_CPDNPD= P(23);           % NPDL1 to PDL1 conversion

K_CSTD=P(24);              % Death rate of Tumor stem

K_CNPDL1D=P(25);           % Death rate of NPDL1

K_CPDL1D=P(26);            % Death rate of PDL1

K_ResD=P(27);              % Degradation of Resources

Delta=P(28);               % Width of the barrier

%% T cell proliferation, interactions, and death

Y_TKM=P(29);               % Carrying capacity for PD1+ T cells

K_TKPD=P(30);              % Proliferation rate of PD1+ T cells 

K_TKNPD=P(31);             % Proliferation rate of PD1- T cells

K_TH=P(32);                % Proliferation rate of helper T cells

K_TREG=P(33);              % Proliferation rate of regulator T cells

K_TEX=P(34);               % Proliferation rate of exhausted T cells

K_TKPDNPD1=P(35);          % Conversion from PD1+ to PD1- T cells

K_TKPDTEX=P(36);           % Conversion from PD1+ T cels to exhausted T cells  

K_THTK=P(37);              % Helper induced proliferation of Killer T cells

K_IL2TK=P(38);             % IL-2 induced proliferation of killer T cells

K_CNPDT=P(39);             % MHC sensing of tumor cells for immune activation

K_CAFTREG=P(40);           % CAF induced proliferation of regulatory T cells

K_TREGTH=P(41);            % TREG induced inhibition of helper T cells  

K_TKPDD=P(42);             % Death rate of killer PD1+ T cells

K_THD=P(43);               % Death rate of helper T cells

K_TREGD=P(44);             % Death rate of regulatory T cells

K_TEXD=P(45);              % Death rate for exhausted T cells

%% CAF proliferation, interactions, and death

Y_FWTM=P(46);              % maximum carrying capacity of Fibroblasts

K_FWT=P(47);               % Proliferation rate of wild type fibroblasts

K_CAF=P(48);               % Proliferation rate of CAF cells

K_OPNCAF=P(49);            % OPN-induced Proliferation rate of CAF cells

K_CTCAF=P(50);             % Tumor cells-induced growth of CAF

K_CTCAFR=P(51);            % Proximity factor for CAF and resistant tumor cells  

K_LIF=P(52);               % LIF fraction available for conversion
 
K_LIFT=P(53);              % Hill like dissociation constant

K_M2CAF=P(54);             % M2-driven proliferation rate for CAFs 

K_FWTCAF=P(55);            % LIF-driven conversion rate from FWT to CAFs

K_CAFFWT=P(56);            % Conversion from CAF to wild type fibroblasts

K_FWTD=P(57);              % Death rate of wild type fibroblasts

K_CAFD=P(58);              % Death rate of CAFs

%% Macrophages

Y_MM= P(59);               % Carrying capacity

K_M1= P(60);               % Growth rate of M1 phase macrophage

K_M2= P(61);               % Growth rate of M2 phase macrophage

K_CANM1= P(62);            % Tumor cells-driven proliferation of macrophages

K_M2M1= P(63);             % Default conversion rate

K_M1M2ICAM1=P(64);         % ICAM1-driven conversion rate

K_ICAM1=P(65);             % Fraction of ICAM1 in proximity with M1 macrophage.

K_CAFM2= P(66);            % Tumor cells-driven proliferation of macrophages

K_M1D= P(67);              % Death rate of M1 phase macrophage 

K_M2D= P(68);              % Death rate of M2 phase macrophage 

%% Cytokines, Chemokines, and Lactate

%% IL-2
K_TIL2=P(69);              % IL-2 secretion by T cells

K_IL2D=P(70);              % IL-2 degradation


%% IFNG
K_TIFNG=P(71);              % IFNG secretion by T cells

K_OPNIFNG=P(72);            % IFNG inhibition by OPN 

K_IFNGD=P(73);              % IFNG degradation rate


%% ICAM1
K_TICAM1=P(74);              % ICAM1 secretion by T cells

K_ICAM1D=P(75);              % ICAM1 degradation rate


%% OPN
K_CANOPN=P(76);              % OPN secretion by Tumor cells

K_CAFOPN=P(77);              % OPN Secretion by CAF

K_IRFOPN=P(78);              % IRF-driven inhibition of OPN

K_OPND=P(79);                % OPN Degradation rate


%% LIF
K_CANLIF=P(80);              % LIF secretion by Tumor cells

K_CAFLIF=P(81);              % LIF Secretion by CAF

K_LIFD=P(82);                % LIF Degradation rate

%% IL-8
K_CANIL8=P(83);              % IL8 secretion by Tumor cells

K_CAFIL8=P(84);              % IL8 Secretion by CAF

K_M2IL8=P(85);               % IL8 secretion by M2

K_IL8D=P(86);                % LIF Degradation rate


%% IRF8
K_M1IRF8=P(87);              % IRF8 secretion by M1

K_IRF8D=P(88);               % IRF8 Degradation rate


%% Lactate
K_M2Lac=P(89);               % Lac secretion by M2

K_CANLac=P(90);              % Lac secretion by Tumor cells

K_LacD=P(91);                % Lac Degradation rate 


%% cell states, cyto and chemokines
CST= y(1);              % T_exposed Tumor stem cells

CSNT= y(2);             % Non T_exposed Tumor stem cells

CNPDL1=y(3);            % Tumor cells exposed to T-cells without pdl1   

CPDL1=y(4);             % Tumor cells exposed to T-cells with pdl1

CRNPDL1=y(5);           % Tumor cells hidden from T-cells without pdl1

CRPDL1=y(6);            % Tumor cells hidden from T-cells with pdl1

Res=y(7);               % Resource concentration

TKPD1= y(8);            % Killer T cells with pd1

TKNPD1= y(9);           % Killer T cells without pd1

TH= y(10);              % Helper T cells

TREG= y(11);            % Regulatory T cells

TEX= y(12);             % Exhausted T cells

FWT= y(13);             % Wild type/ myofibroblasts

CAF= y(14);             % Cancer associated fibroblasts

MACM1= y(15);           % M1 phase macrophage

MACM2= y(16);           % M2 phase macrophage

IL2= y(17);             % Interleukin 2

LIF= y(18);             % Leukemia inhibitory factor

IFNG= y(19);            % Interferon gamma

IL8= y(20);             % Interleukin 8

LAC= y(21);             % Lactate

ICAM= y(22);            % Inter-cellular adhesion molecule

OPN= y(23);             % Osteopontin    

IRF8= y(24);            % Interferon regulatory factor 8

%% Hyper parameters
I=tanh(alpha*K_CAFB*CAF);

alpha_T=exp(-Delta*alpha^2*CAF^2);

%% Population of cancer cell states

% T-exposed tumor stem cell

F_ResCST=Res*K_RCST*CST*(1-CST/((1-I)*Y_CST+1))*1/(alpha_Comp*(CPDL1+CNPDL1)+1);          % Resource-driven proliferation of tumor stem cells

F_TxCST=(1+K_TXC*TEX/(TEX+1));                                                          % Exhausted T-cells-driven proliferation modulator of tumor stem cells 

F_CAFCST=(1+K_CAFC*CAF/(CAF+1));                                                        % OPN-driven proliferation modulator of tumor stem cells

F_CSTCNPDL1= K_CSTCNPDL1*(IL8/(IL8+1))*CST;                                             % Conversion from stem to non-PDL1, T-exposed tumor cells

F_CSTCPDL1= K_CSTCPDL1*(IL8/(IL8+1))*CST;                                               % Conversion from stem to PDL1, T-exposed tumor cells

F_TKCST= K_TKC*(IFNG/(IFNG+1))*CST*(TKPD1+TKNPD1)*1/(K_Lac*LAC+1);

F_DCST= K_CSTD*(1+IFNG/(IFNG+1)+MACM1/(MACM1+1))*CST;                                   % Death of stem cells

dydt_CST=F_ResCST*F_CAFCST*F_TxCST-F_TKCST-(F_CSTCNPDL1+F_CSTCPDL1)-F_DCST;

% Non T_exposed Tumor Stem cell

F_ResCSNT=Res*K_RCSNT*CSNT*(1-CSNT/(I*Y_CST+1))*1/(alpha_Comp*(CRPDL1+CRNPDL1)+1);        % Resource-driven proliferation of tumor stem cells

F_CAFCSNT=(1+K_CAFCR*CAF/(1+CAF));                                                      % CAF-driven proliferation modulator of tumor stem cells

F_TxCSNT=(1+K_TXC*TEX*alpha_T/(TEX+1));                                                 % Exhausted T-cells-driven proliferation modulator of tumor stem cells 

F_TKCSNT= K_TKC*(IFNG/(IFNG+1))*CSNT*alpha_T*(TKPD1+TKNPD1)*1/(1+K_Lac*LAC);

F_CSNTCRNPDL1=K_CSTCNPDL1*IL8/(IL8+1)*CSNT;                                             % Conversion from stem to non-PDL1, non T-exposed tumor cells

F_CSNTCRPDL1= K_CSTCPDL1*IL8/(IL8+1)*CSNT;                                              % Conversion from stem to PDL1, non T-exposed tumor cells

F_DCSNT= K_CSTD*(1+IFNG/(IFNG+1)+MACM1/(MACM1+1))*CSNT;                                 % Death of stem cells

dydt_CSNT=F_ResCSNT*F_CAFCSNT*F_TxCSNT-F_TKCSNT-(F_CSNTCRNPDL1+F_CSNTCRPDL1)-F_DCSNT;


% PDL1-, T-exposed tumor cells

F_ResCNPDL1= Res*K_RCNPDL1*CNPDL1*(1-CNPDL1/((1-I)*Y_CPDL1M+1))*1/(alpha_Comp*(CST+CPDL1)+1);  % Resource-driven proliferation of tumor stem cells

F_CAFCNPDL1= (1+K_CAFC*CAF/(1+CAF));                                                         % CAF-based proliferation modulator of PDL1-ve tumor cells
 
F_TKCNPDL1=  K_TKC*IFNG/(IFNG+1)*(TKPD1+TKNPD1)*CNPDL1*1/(1+LAC*K_Lac);                      % Killer T-cell driven death of PDL1-ve tumor cells

F_CNPDL1CPDL1= K_CPDNPD*IFNG/(IFNG+1)*CNPDL1;                                                % IFNG-induced conversion from PDL1-ve to +ve tumor cells

F_DCNPDL1= K_CNPDL1D*CNPDL1;                                                                 % Death of Non-PDL1, T-exposed tumor cells

dydt_CNPDL1= F_ResCNPDL1*F_CAFCNPDL1*F_TxCST+F_CSTCNPDL1-F_TKCNPDL1-F_CNPDL1CPDL1-F_DCNPDL1;


% PDL1+, T-exposed tumor cells

F_ResCPDL1= Res*K_RCPDL1*CPDL1*(1-CPDL1/(Y_CPDL1M*(1-I)+1))*1/(alpha_Comp*(CST+CNPDL1)+1);     % Resource-driven proliferation of PDL1 +ve cells

F_CAFCPDL1=(1+K_CAFC*CAF/(1 + CAF));                                                         % CAF-based proliferation modulator of PDL1+ve tumor cells

F_TKNPDCPDL1=K_TKC*IFNG/(IFNG+1)*CPDL1*TKNPD1*1/(K_Lac*LAC+1);                               % Killer, PD1-ve T-cell driven death of PDL1+ve tumor cells

F_DCPDL1= K_CPDL1D*CPDL1;                                                                    % Death of PDL1, T-exposed tumor cells   

dydt_CPDL1=F_ResCPDL1*F_CAFCPDL1*F_TxCST+F_CNPDL1CPDL1+F_CSTCPDL1-F_TKNPDCPDL1-F_DCPDL1;


% PDL1-, non-T-exposed tumor cells 

F_ResCRNPDL1=Res*K_RCRNPDL1*CRNPDL1*(1-CRNPDL1/(I*Y_CPDL1M+1))*1/(alpha_Comp*(CSNT+CRPDL1)+1); % Resource-driven proliferation of Non-PDL1, Non T-exposed tumor cells

F_CAFCRNPDL1=(1+K_CAFCR*CAF/(1 + CAF));                                                      % CAF-based proliferation modulator of Non-PDL1, Non T-exposed PDL1+ve tumor cells

F_TKCRNPDL1= K_TKC*IFNG/(IFNG+1)*alpha_T*(TKPD1+TKNPD1)*CRNPDL1*1/(K_Lac*LAC+1);             % Killer T-cell driven death of PDL1-ve tumor cells

F_CRNPDL1CRPDL1=  K_CPDNPD*IFNG/(IFNG+1)*CRNPDL1;                                            % ING-induced conversion from PDL1-ve to +ve tumor cells

F_DCRNPDL1= K_CNPDL1D*CRNPDL1;                                                               % Death rate of non T-exposed, non-PDL1 tumor cells


dydt_CRNPDL1= F_ResCRNPDL1*F_CAFCRNPDL1*F_TxCSNT+F_CSNTCRNPDL1-F_TKCRNPDL1-F_CRNPDL1CRPDL1-F_DCRNPDL1;


% PDL1+, non T-exposed tumor cells

F_ResCRPDL1= Res*K_RCRPDL1*CRPDL1*(1-CRPDL1/(I*Y_CPDL1M+1))*1/(alpha_Comp*(CSNT+CNPDL1)+1);     % Resource-driven proliferation of Non-PDL1, Non T-exposed tumor cells

F_CAFCRPDL1=(1+K_CAFCR*CAF/(1 + CAF));                                                        % CAF-based proliferation modulator of PDL1+ve tumor cells

F_TKNPDCRPDL1= K_TKC*IFNG/(IFNG+1)*TKNPD1*CRPDL1*1/(LAC*K_Lac+1)*alpha_T;                     % Killer, PD1-ve T-cell driven death of PDL1+ve tumor cells

F_DCRPDL1= K_CPDL1D*CRPDL1;                                                                   % Death of PDL1, T-exposed tumor cells   


dydt_CRPDL1=F_ResCRPDL1*F_CAFCRPDL1*F_TxCSNT+F_CRNPDL1CRPDL1+F_CSNTCRPDL1-F_TKNPDCRPDL1-F_DCRPDL1;


%% Resource concentration

F_ProRes=K_RIn*(Y_RM-Res);

F_DRes= K_ResD*Res;

dydt_Res=F_ProRes-F_DRes;

%% T-cell population

% Killer PD1+ T cells

F_ProTKPD1=K_TKPD*TKPD1*(1-TKPD1/(Y_TKM-TEX-TKNPD1+1));                                   % Proliferation of killer, PD1+ve T-cells

F_THTKPD1=(1+K_THTK*(TH*MACM1)/(TH*MACM1+1));                                           % T-Helper-driven proliferation modulator of Killer, PD1+ve T cells

F_IL2TKPD1=(1+K_IL2TK*IL2/(1+IL2));                                                     % IL2-driven proliferation modulator of Killer, PD1+ve T cells

F_TKPD1TKNPD1= K_TKPDNPD1*u*TKPD1;                                                      % ICI-based conversion from PD1+ve to -ve killer T cells

F_TKPD1TEX=K_TKPDTEX*TKPD1*(CPDL1*MACM2)/(CPDL1*MACM2+1);                               % PDL1+ve and M2 driven conversion of PD1+ve to exhausted T cells  

F_DTKPD1= K_TKPDD*TKPD1;                                                                % Death rate of PD1+, killer T cells

dydt_TKPD1 = F_ProTKPD1*F_THTKPD1*F_IL2TKPD1 - F_TKPD1TKNPD1 - F_TKPD1TEX - F_DTKPD1; 


% Killer PD1- T cells
F_ProTKNPD1=K_TKNPD*TKNPD1*u*(1-TKNPD1/(Y_TKM-TEX-TKPD1+1));                              % Proliferation of killer, PD1+ve T-cells

F_THTKNPD1=(1+u*K_THTK*TH/(TH+1));                                                      % T-Helper-driven proliferation modulator of Killer, PD1+ve T cells

F_IL2TKNPD1=(1+K_IL2TK*IL2/(1+IL2));                                                    % IL2-driven proliferation modulator of Killer, PD1+ve T cells

F_DTKNPD1= K_TKPDD*TKNPD1;

%dydt_TKNPD1=F_TKPD1TKNPD1-F_DTKNPD1;

dydt_TKNPD1=F_ProTKNPD1*F_THTKNPD1*F_IL2TKNPD1+F_TKPD1TKNPD1-F_DTKNPD1;


% Helper T cells

F_ProTH=K_TH*TH*(1-TH/Y_TKM);                                                           % Proliferation of helper T-cells

F_CNPDL1TH=(1+K_CNPDT*(CNPDL1+u*CPDL1)/(CNPDL1+u*CPDL1+1));                             % PDL1-ve tumor cell-driven proliferation modulator of helper T cells

F_TREGTH=1/(1+K_TREGTH*TREG);                                                           % Regulatory T cell-driven regulation of helper T cell population

F_DTH=K_THD*TH;                                                                         % Death rate of Helper T cells

dydt_TH=F_ProTH*F_CNPDL1TH*F_TREGTH-F_DTH;


% Regulation T cells

F_ProTREG=K_TREG*TREG*(1-TREG/Y_TKM);                                                   % Proliferation of regulatory T cells

F_CAFTREG=(1+K_CAFTREG*CAF/(1+CAF));                                                    % CAF-driven proliferation of regulatory T cells                  

F_DTREG=K_TREGD*TREG;                                                                   % Death rate of Regulatory T cells

dydt_TREG=F_ProTREG*F_CAFTREG-F_DTREG;



% Exhausted T cells
F_ProTEX=K_TEX*TEX*(1-TEX/(Y_TKM-TKPD1-TKNPD1+1));

F_DTEX=K_TEXD*TEX;                                                                      % Death rate of exhausted t cells

dydt_TEX=F_ProTEX+F_TKPD1TEX-F_DTEX;



%% Fibroblasts population

% Wild type

F_ProFWT=K_FWT*FWT*(1-FWT/(Y_FWTM-CAF+1));                                                % Proliferation of Wild-type fibroblasts

F_CAFFWT=K_CAFFWT*CAF;                                                                  % Conversion from CAF to wild type fibroblasts

F_FWTCAF=K_FWTCAF*(K_LIF*LIF)^2/((K_LIF*LIF)^2+K_LIFT)*FWT;                             % LIF-driven conversion from wild type fibroblasts to CAF

F_DFWT=K_FWTD*FWT;                                                                      % Death rate of wild type fibroblast

dydt_FWT=F_ProFWT-F_FWTCAF+F_CAFFWT-F_DFWT;


% Invasive

F_ProCAF=K_CAF*CAF*(1-CAF/(Y_FWTM-FWT+1));                                               % Proliferation of CAFs

F_OPNCAF=(1+K_OPNCAF*OPN/(1+OPN));                                                     % OPN-driven proliferation modulator of CAFs

F_M2CAF= (1+K_M2CAF*MACM2/(1+MACM2));                                                  % OPN-driven proliferation modulator of CAFs

F_CANCAF=(1+K_CTCAF*((CST+CNPDL1+CPDL1)/(1+CST+CNPDL1+CPDL1)+K_CTCAFR*(CSNT+CRNPDL1+CRPDL1)/(1+CSNT+CRNPDL1+CRPDL1)));

F_DCAF=K_CAFD*CAF;                                                                     % Death rate of cancer associated fibroblasts

dydt_CAF=F_ProCAF*F_OPNCAF*F_OPNCAF*F_M2CAF*F_CANCAF+F_FWTCAF-F_CAFFWT-F_DCAF;

%% Macrophages

% M1 phase

F_ProMACM1= K_M1*MACM1*(1-MACM1/(Y_MM-MACM2+1));                                         % Proliferation of M1 macrophage

F_CANM1= (1+K_CANM1)*(CNPDL1+u*CPDL1+CST+CSNT+u*CRPDL1+CRNPDL1)/(CNPDL1+u*CPDL1+CST+CSNT+u*CRPDL1+CRNPDL1+1);

F_M1M2= K_M1M2ICAM1*MACM1*(K_ICAM1*ICAM)^2/((K_ICAM1*ICAM)^2+1);                       % Lactate-driven Conversion from M1 to M2 macrophage

F_M2M1= K_M2M1*MACM2;                                                                  % Conversion from M2 to M1 macrophage

F_DMACM1= K_M1D*MACM1;                                                                 % Death rate of M1 macrophage

dydt_MACM1=F_ProMACM1*F_CANM1+F_M2M1-F_M1M2-F_DMACM1;

% M2 phase

F_ProMACM2=K_M2*MACM2*(1-MACM2/(Y_MM-MACM1+1));                                          % Proliferation of M2 macrophage

F_CAFMACM2=(1+K_CAFM2*CAF/(CAF+1));

F_DMACM2= K_M2D*MACM2;                                                                 % Death rate of M2 macrophage

dydt_MACM2= F_ProMACM2*F_CAFMACM2+F_M1M2-F_M2M1-F_DMACM2;

%% Cytokines concentration

% Interleukin 2

F_TKIL2=K_TIL2*(TKPD1+TKNPD1);                                                         % Killer T-cells secretion of IL-2

F_CIL2=C_IL2;                                                                          % IL-2 injection

F_DIL2=K_IL2D*IL2;                                                                     % Degradation rate of IL2

dydt_IL2=F_CIL2+F_TKIL2-F_DIL2;

% Lukemia Inhibiting Factor (IL6)

F_CANLIF=(CST+CSNT+CNPDL1+CPDL1+CRNPDL1+CRPDL1)*K_CANLIF;                              % CAF secretion of LIF

F_CAFLIF=(CAF)*K_CAFLIF;                                                               % IFNG-driven inhibition of LIF

F_DLIF=K_LIFD*LIF;                                                                     % Death rate of LIF

F_LIFKnock=C_LIF*LIF;

dydt_LIF=F_CAFLIF+F_CANLIF-F_DLIF-F_LIFKnock;


% Interferon Gamma

F_TKIFN=K_TIFNG*(TKPD1+TKNPD1);                                                        % Killer T-cell secretion of IFNG

F_OPNIFN=1/(OPN*K_OPNIFNG+1);                                                          % OPN-driven inhibition of IFNG

F_DIFNG=K_IFNGD*IFNG;                                                                  % Death rate of IFNG

dydt_IFNG=F_TKIFN*F_OPNIFN-F_DIFNG;

% Interleukin 8

F_CAFIL8= K_CAFIL8*CAF;                                                                % CAF secretion of IL-8

F_CANIL8= K_CANIL8*(CST+CSNT+CPDL1+CNPDL1+CRPDL1+CRNPDL1);

F_M2IL8=  K_M2IL8*MACM2;

F_DIL8=K_IL8D*IL8;                                                                     % Death rate of IL-8

F_IL8Knock= C_IL8*IL8;

dydt_IL8=F_CAFIL8+F_CANIL8+F_M2IL8-F_IL8Knock-F_DIL8;

% Lactate

F_CANLac=  (CST+CSNT+CPDL1+CNPDL1+CRPDL1+CRNPDL1)*K_CANLac;                           % Tumor cells secretion of Lactate

F_MACLac=  K_M2Lac*MACM2;                                                             % M2 macrophage secretion of Lactate

F_LacKnock= C_Lac*LAC;

F_DLAC= K_LacD*LAC;                                                                   % Death rate of LAC


dydt_LAC=F_CANLac+F_MACLac-F_LacKnock-F_DLAC;



% Intercellular Adhesion Molecule

F_TKICAM=K_TICAM1*(TKPD1+TKNPD1);                                                     % T cell-based secretion of ICMA

F_DICAM= K_ICAM1D*ICAM;                                                               % Death rate of ICMA

dydt_ICAM= F_TKICAM-F_DICAM;



% Osteopontin

F_CAFOPN=K_CAFOPN*CAF;                                                                % CAF secretion of OPN

F_CANOPN=K_CANOPN*(CST+CSNT+CPDL1+CNPDL1+CRPDL1+CRNPDL1);                             % Tumor cells secretion of OPN

F_IRF8OPN=1/(K_IRFOPN*IRF8+1);                                                        % IRF8-based inhibition of OPN

F_DOPN=K_OPND*OPN;                                                                    % Degradation rate of OPN

F_OPNKnock=C_OPN*OPN;

dydt_OPN=(F_CAFOPN+F_CANOPN)*F_IRF8OPN-F_DOPN-F_OPNKnock;



% Interferon Regulatory Factor 8

F_MACM1IRF8=K_M1IRF8*MACM1;                                                            % M1 macrophage secretion of IRF8

F_DIRF8=K_IRF8D*IRF8;                                                                  % Degradation rate of IRF8

dydt_IRF8= F_MACM1IRF8-F_DIRF8;


dydt=[dydt_CST;dydt_CSNT;dydt_CNPDL1;dydt_CPDL1;dydt_CRNPDL1;dydt_CRPDL1;dydt_Res;dydt_TKPD1;dydt_TKNPD1;dydt_TH;dydt_TREG;dydt_TEX;dydt_FWT;dydt_CAF;dydt_MACM1;dydt_MACM2;dydt_IL2;dydt_LIF;dydt_IFNG;dydt_IL8;dydt_LAC;dydt_ICAM;dydt_OPN;dydt_IRF8];

end

