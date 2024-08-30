function P=HNSCC_parameters_OPN(alpha)

%% Load the parameters

%% Resource to cancer cells
K_RCST=100;
P(1)=K_RCST;               %  per unit of resource per unit of time. Resource to T-exposed tumor stem cell 

K_RCSNT=100;
P(2)=K_RCSNT;              %  per unit of resource per unit of time. Resource to Non-T-exposed tumor stem cell

K_RCNPDL1=120;            
P(3)=K_RCNPDL1;            %  per unit of resource per unit of time. Resource to T-exposed PDL1- tumor cell 

K_RCRNPDL1=120;           
P(4)=K_RCRNPDL1;           %  per unit of resource per unit of time. Resource to Non-T-exposed PDL1- tumor cell

K_RCPDL1=10;
P(5)=K_RCPDL1;             %  per unit of resource per unit of time. Resource to T-exposed PDL1+ tumor cell

K_RCRPDL1=10;
P(6)=K_RCRPDL1;            %  per unit of resource per unit of time. Resource to Non-T-exposed PDL1+ tumor cell

K_RIn=1500;
P(7)=K_RIn;                %  Resource per unit time. Default resource intake rate

alpha_Comp=0.0005;               % Resource consumption rate
P(8)=alpha_Comp;
%% Cancer cell proliferation, interactions, and death
Y_CST=10^4;
P(9)=Y_CST;            % Carrying capacity for PD1- T cells

Y_CPDL1M=10^4;
P(10)=Y_CPDL1M;            % Carrying capacity for PD1+ T cells

Y_RM=2;
P(11)=Y_RM;                % Resources. Maximum resource withholding capacity.

K_TXC=50;
P(12)=K_TXC;               % Exhausted T cells-driven growth of tumor cells.

K_FWTC=0.005;
P(13)=K_FWTC;              % Fraction of wild type fibroblasts in proximity with the tumor cells.

K_CAFC=80;
P(14)=K_CAFC;              % CAF-driven growth of the tumor cells.

K_CAFCR=2.5;
P(15)=K_CAFCR;             % CAF-driven growth of the tumor non-T exposed cells.
            
K_TKC=1500;
P(16)=K_TKC;               % T cells-driven apoptosis of tumor cells.

K_Lac=0.05;
P(17)=K_Lac;               % Inhibition rate of lactate

K_TKCIFNG=1;
P(18)=K_TKCIFNG;           % IFNG induced T cells-driven apoptosis of tumor cells.

P(19)=alpha;               % CAF barrier

K_CAFB=0.001;
P(20)=K_CAFB;              % Barrier formation rate

K_CSTCNPDL1=8;
P(21)=K_CSTCNPDL1;         % Stem to PDL1- 

K_CSTCPDL1=8;
P(22)=K_CSTCPDL1;          % Stem to PDL1+ 

K_CPDNPD=8;
P(23)=K_CPDNPD;            % NPDL1 to PDL1 conversion

K_CSTD=0.0005;
P(24)=K_CSTD;              % Death rate of Tumor stem

K_CNPDL1D=8;
P(25)=K_CNPDL1D;           % Death rate of NPDL1

K_CPDL1D=8;
P(26)=K_CPDL1D;            % Death rate of PDL1

K_ResD=8;
P(27)=K_ResD;              % Degradation of Resources

Delta=10^-3;
P(28)=Delta;               % Width of the barrier

%% T cell proliferation, interactions, and death

Y_TKM=5000;
P(29)=Y_TKM;               % Carrying capacity for PD1+ T cells

K_TKPD=45;
P(30)=K_TKPD;              % Proliferation rate of PD1+ T cells 

K_TKNPD=60;
P(31)=K_TKNPD;             % Proliferation rate of PD1- T cells

K_TH=25;
P(32)=K_TH;                % Proliferation rate of helper T cells

K_TREG=10;
P(33)=K_TREG;              % Proliferation rate of regulator T cells

K_TEX=15;
P(34)=K_TEX;               % Proliferation rate of Exhausted T cells

K_TKPDNPD1=150000;
P(35)=K_TKPDNPD1;          % Conversion from PD1+ to PD1- T cells

K_TKPDTEX=15;
P(36)=K_TKPDTEX;           % Conversion from PD1+ T cels to exhausted T cells  

K_THTK=40;
P(37)=K_THTK;              % Helper induced proliferation of Killer T cells

K_IL2TK=20;
P(38)=K_IL2TK;             % IL-2 induced proliferation of killer T cells

K_CNPDT=90;
P(39)=K_CNPDT;             % MHC sensing of tumor cells for immune activation

K_CAFTREG=15;
P(40)=K_CAFTREG;           % CAF induced proliferation of regulatory T cells

K_TREGTH=0.001;
P(41)=K_TREGTH;            % TREG induced inhibition of helper T cells  

K_TKPDD=10;
P(42)=K_TKPDD;             % Death rate of killer PD1+ T cells

K_THD=7;
P(43)=K_THD;               % Death rate of helper T cells

K_TREGD=9;
P(44)=K_TREGD;             % Death rate of regulatory T cells

K_TEXD=9;
P(45)=K_TEXD;              % Death rate for exhausted T cells

%% CAF proliferation, interactions, and death

Y_FWTM=5000;
P(46)=Y_FWTM;              % maximum carrying capacity of Fibroblasts

K_FWT=50;
P(47)=K_FWT;               % Proliferation rate of wild type fibroblasts

K_CAF=1;
P(48)=K_CAF;               % Proliferation rate of CAF cells

K_OPNCAF=50;
P(49)=K_OPNCAF;            % OPN-induced Proliferation rate of CAF cells

K_CTCAF=0.05;
P(50)=K_CTCAF;             % Tumor cells-induced growth of CAF

K_CTCAFR=1;
P(51)=K_CTCAFR;            % Proximity factor for CAF and resistant tumor cells  

K_LIF=0.0005;
P(52)=K_LIF;               % LIF fraction available for conversion

K_LIFT=1;              
P(53)=K_LIFT;              % Hill like dissociation constant

K_M2CAF=0.05;
P(54)=K_M2CAF;             % M2-driven proliferation rate for CAFs 

K_FWTCAF=50;
P(55)=K_FWTCAF;            % LIF-driven conversion rate from FWT to CAFs

K_CAFFWT=9;
P(56)=K_CAFFWT;            % Conversion from CAF to wild type fibroblasts

K_FWTD=10;
P(57)=K_FWTD;              % Death rate of wild type fibroblasts

K_CAFD=15;
P(58)=K_CAFD;              % Death rate of CAFs

%% Macrophage proliferation, interactions, and death

Y_MM=5000;
P(59)=Y_MM;                % Carrying capacity

K_M1=15;
P(60)=K_M1;                % Growth rate of M1 phase macrophage

K_M2=5;
P(61)=K_M2;                % Growth rate of M2 phase macrophage

K_CANM1=10;
P(62)=K_CANM1;             % Tumor cells-driven proliferation of macrophages

K_M2M1=3;
P(63)=K_M2M1;              % Default conversion rate

K_M1M2ICAM1=30;
P(64)=K_M1M2ICAM1;         % ICAM1-driven conversion rate

K_ICAM1=5*10^-4;
P(65)=K_ICAM1;             % Fraction of ICAM1 in proximity with M1 macrophage.


K_CAFM2=120;
P(66)=K_CAFM2;             % Tumor cells-driven proliferation of macrophages

K_M1D=10;
P(67)=K_M1D;               % Death rate of M1 phase macrophage 

K_M2D=10;
P(68)=K_M2D;               % Death rate of M2 phase macrophage 

%% Cytokines, Chemokines, and Lactate
%% IL-2
K_TIL2=5;
P(69)=K_TIL2;               % IL-2 secretion by T cells

K_IL2D=6;
P(70)=K_IL2D;               % IL-2 degradation


%% IFNG
K_TIFNG=20;
P(71)=K_TIFNG;              % IFNG secretion by T cells

K_OPNIFNG=0.01;
P(72)=K_OPNIFNG;            % IFNG inhibition by OPN 

K_IFNGD=5;
P(73)=K_IFNGD;              % IFNG degradation rate


%% ICAM1
K_TICAM1=5;
P(74)=K_TICAM1;             % ICAM1 secretion by T cells

K_ICAM1D=5;
P(75)=K_ICAM1D;             % ICAM1 degradation rate


%% OPN
K_CANOPN=3;
P(76)=K_CANOPN;             % OPN secretion by Tumor cells

K_CAFOPN=5;
P(77)=K_CAFOPN;             % OPN Secretion by CAF

K_IRFOPN=0.05;
P(78)=K_IRFOPN;             % IRF-driven inhibition of OPN

K_OPND=4;
P(79)=K_OPND;               % OPN Degradation rate


%% LIF
K_CANLIF=6;
P(80)=K_CANLIF;             % LIF secretion by Tumor cells

K_CAFLIF=0.05;
P(81)=K_CAFLIF;             % LIF Secretion by CAF

K_LIFD=8;
P(82)=K_LIFD;               % LIF Degradation rate

%% IL-8
K_CANIL8=2;
P(83)=K_CANIL8;             % IL8 secretion by Tumor cells

K_CAFIL8=2;
P(84)=K_CAFIL8;             % IL8 Secretion by CAF

K_M2IL8=15;
P(85)=K_M2IL8;              % IL8 secretion by M2

K_IL8D=5;
P(86)=K_IL8D;               % LIF Degradation rate


%% IRF8
K_M1IRF8=2;
P(87)=K_M1IRF8;             % IRF8 secretion by M1

K_IRF8D=2;
P(88)=K_IRF8D;              % IRF8 Degradation rate


%% Lactate
K_M2Lac=2;
P(89)=K_M2Lac;              % Lac secretion by M2

K_CANLac=0.2;
P(90)=K_CANLac;             % Lac secretion by Tumor cells

K_LacD=4;
P(91)=K_LacD;               % Lac Degradation rate 

end