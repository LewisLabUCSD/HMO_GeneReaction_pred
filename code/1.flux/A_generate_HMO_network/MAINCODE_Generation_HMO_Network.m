%% STEP1 - Generation of the initial reaction for Human Milk
%% Oligosaccharide synthesis
    display 'Step1 - Initialize HMO Network'
    BasicNetwork = createModel({'HMO_Init'},{'HMO_Init'},{' ->	Ab4G;HMO[c]'},false,1,1,{'HMO Synthesis'});

    BasicNetwork.c(1) = 1;
    BasicNetwork.genes = {'N/A'};
    BasicNetwork.geneLoc = {'N/A'};
    BasicNetwork.rxnGeneMat = true(length(BasicNetwork.rxns),1);
    BasicNetwork.description = 'Basic HMO Network';

    %pause;

%% STEP 2 - Systematic generation of all the oligosaccharide structures presenting
%a defined complexity level (fixed as the maximum number of sugar units of measured HMOs) based on reaction laws  
    display 'Step2 - Create generic HMO Network'    
    K=7;
    HMO_Network = CreateHMONetwork(BasicNetwork,K)
    
    %pause;
    
%% STEP 3 - Definition of HMO candidates structures (5 HMOs do not have
%% defined structures)
    display 'Definition of candidate HMO structures'
    %Known structures
    HMO_0= 'Ab4G;HMO';
    HMO_01 = 'Fa2Ab4G;HMO'; %2FL
    HMO_02 = 'Ab4(Fa3)G;HMO';%3FL
    HMO_03 = 'NNa3Ab4G;HMO'; %3SL
    HMO_04 = 'Ab3GNb3Ab4G;HMO';%LNT
    HMO_05 = 'Ab4GNb3Ab4G;HMO';%LNnT
    HMO_06 = 'Fa2Ab3GNb3Ab4G;HMO';%LNFPI 
    HMO_07 = 'Ab3(Fa4)GNb3Ab4G;HMO';% LNFPII
    HMO_08 = 'Ab4(Fa3)GNb3Ab4G;HMO';% LNFPIII
    HMO_09 = 'Ab3(NNa6)GNb3Ab4G;HMO';% LSTb  
    HMO_10 = 'NNa6Ab4GNb3Ab4G;HMO';% LSTc 
    HMO_11 = 'NNa3Ab3(NNa6)GNb3Ab4G;HMO';%DSLNT 
    
    %Unknown structures
    HMO_12 = 'Fa2Ab3(Fa4)GNb3Ab4G;HMO';% DFLNT1 
    HMO_13 = 'Ab3(Fa4)GNb3Ab4(Fa3)G;HMO';% DFLNT2 
    HMO_14 = 'Ab4(Fa3)GNb3Ab4(Fa3)G;HMO';% DFLNT3 
    
    HMO_15 = 'Ab3GNb3Ab4(Fa3)GNb3Ab4G;HMO';% FLNH1 	
    HMO_16 = 'Ab3(Fa4)GNb3Ab4GNb3Ab4G;HMO';% FLNH2 
    HMO_17 = 'Fa2Ab3GNb3(Ab4GNb6)Ab4G;HMO';% FLNH3 	
    HMO_18 = 'Ab3GNb3(Ab4(Fa3)GNb6)Ab4G;HMO';% FLNH4 
   
    HMO_19 = 'Ab3(Fa4)GNb3(Ab4GNb6)Ab4G;HMO';% FLNH5 
    HMO_20 = 'Ab4GNb3Ab4(Fa3)GNb3Ab4G;HMO';% FLNH6 	
    HMO_21 = 'Fa2Ab3GNb3Ab4GNb3Ab4G;HMO';% FLNH7 
    
    
    HMO_22 = 'Fa2Ab3GNb3(Ab4(Fa3)GNb6)Ab4G;HMO';% DFLNH1  
    HMO_23 = 'Ab3(Fa4)GNb3(Ab4(Fa3)GNb6)Ab4G;HMO';% DFLNH2  
    HMO_24 = 'Ab3(Fa4)GNb3Ab4(Fa3)GNb3Ab4G;HMO';% DFLNH3  	
    HMO_25 = 'Fa2Ab3GNb3Ab4(Fa3)GNb3Ab4G;HMO';% DFLNH4  
    HMO_26 = 'Ab4(Fa3)GNb3Ab4(Fa3)GNb3Ab4G;HMO';% DFLNH5  	
    HMO_27 = 'Ab4(Fa3)GNb3(Ab4(Fa3)GNb6)Ab4G;HMO';% DFLNH6  
    HMO_28 = 'Fa2Ab3(Fa4)GNb3(Ab4GNb6)Ab4G;HMO';% DFLNH7 
    
    
    HMO_29 = 'NNa3Ab3(NNa6)GNb3(Fa2Ab4GNb6)Ab4G;HMO';% FDSLNH1	
    HMO_30 = 'NNa3Ab3(NNa6)GNb3(Ab4(Fa3)GNb6)Ab4G;HMO';% FDSLNH2	
    HMO_31 = 'NNa3Ab3(Fa4)GNb3(NNa6Ab4GNb6)Ab4G;HMO';% FDSLNH3	
    HMO_32 = 'NNa6Ab4(Fa3)GNb3(NNa6Ab4GNb6)Ab4G;HMO';% FDSLNH4		
    HMO_33 = 'NNa6Ab4(Fa3)GNb3(NNa3Ab4GNb6)Ab4G;HMO';% FDSLNH5		
    HMO_34 = 'NNa3Ab4(Fa3)GNb3(NNa6Ab4GNb6)Ab4G;HMO';% FDSLNH6	
    HMO_35 = 'NNa3Ab4(Fa3)GNb3(NNa3Ab4GNb6)Ab4G;HMO';% FDSLNH7	
    
    HMO_36 = 'NNa3Ab3GNb3(NNa6Ab4GNb6)Ab4G;HMO';% DSLNH1  
    HMO_37 = 'NNa3Ab3(NNa6)GNb3(Ab4GNb6)Ab4G;HMO';% DSLNH2  
    HMO_38 = 'NNa6Ab4GNb3(NNa6Ab4GNb6)Ab4G;HMO';% DSLNH3 
    
    HMO_39= 'DFLNT;HMO';
    HMO_40= 'FLNH;HMO';
    HMO_41= 'DFLNH;HMO';
    HMO_42= 'FDSLNH;HMO';
    HMO_43= 'DSLNH;HMO';

    WT_Profile.List = {HMO_01,HMO_02,HMO_03,HMO_04,HMO_05,HMO_06,HMO_07,...
    HMO_08,HMO_09,HMO_10,HMO_11,HMO_39,HMO_40,HMO_41,HMO_42,HMO_43};
    WT_Profile.Fractions = [1/16.*ones(1,16)];

    %pause;


%% STEP 4 - Addition of 28 reactions to link the 23 candidates structures
%% associated with the 5 measured HMOs having no defined structures
    display 'Addition of reactions to link unknown HMO candidate structures'
    
    [HMO_Network,RxnCheck] = addReaction(HMO_Network,'DFLNT1_rxn',{'Fa2Ab3(Fa4)GNb3Ab4G;HMO[c]','DFLNT;HMO[c]'},[-1 1],false,0,1000,0,'HMO Synthesis');
    [HMO_Network,RxnCheck] = addReaction(HMO_Network,'DFLNT2',{'Ab3(Fa4)GNb3Ab4(Fa3)G;HMO[c]','DFLNT;HMO[c]'},[-1 1],false,0,1000,0,'HMO Synthesis');
    [HMO_Network,RxnCheck] = addReaction(HMO_Network,'DFLNT3',{'Ab4(Fa3)GNb3Ab4(Fa3)G;HMO[c]','DFLNT;HMO[c]'},[-1 1],false,0,1000,0,'HMO Synthesis');

    [HMO_Network,RxnCheck] = addReaction(HMO_Network,'HMOSecr_DFLNT',['DFLNT;HMO[c] ->'],[],false,0,1000,0,'HMO Secretion');

    [HMO_Network,RxnCheck] = addReaction(HMO_Network,'FLNH1',{'Ab3GNb3Ab4(Fa3)GNb3Ab4G;HMO[c]','FLNH;HMO[c]'},[-1 1],false,0,1000,0,'HMO Synthesis');
    [HMO_Network,RxnCheck] = addReaction(HMO_Network,'FLNH2',{'Ab3(Fa4)GNb3Ab4GNb3Ab4G;HMO[c]','FLNH;HMO[c]'},[-1 1],false,0,1000,0,'HMO Synthesis');
    [HMO_Network,RxnCheck] = addReaction(HMO_Network,'FLNH3',{'Fa2Ab3GNb3(Ab4GNb6)Ab4G;HMO[c]','FLNH;HMO[c]'},[-1 1],false,0,1000,0,'HMO Synthesis');
    [HMO_Network,RxnCheck] = addReaction(HMO_Network,'FLNH4',{'Ab3GNb3(Ab4(Fa3)GNb6)Ab4G;HMO[c]','FLNH;HMO[c]'},[-1 1],false,0,1000,0,'HMO Synthesis');
    
    [HMO_Network,RxnCheck] = addReaction(HMO_Network,'FLNH5',{'Ab3(Fa4)GNb3(Ab4GNb6)Ab4G;HMO[c]','FLNH;HMO[c]'},[-1 1],false,0,1000,0,'HMO Synthesis');
    [HMO_Network,RxnCheck] = addReaction(HMO_Network,'FLNH6',{'Ab4GNb3Ab4(Fa3)GNb3Ab4G;HMO[c]','FLNH;HMO[c]'},[-1 1],false,0,1000,0,'HMO Synthesis');
    [HMO_Network,RxnCheck] = addReaction(HMO_Network,'FLNH7',{'Fa2Ab3GNb3Ab4GNb3Ab4G;HMO[c]','FLNH;HMO[c]'},[-1 1],false,0,1000,0,'HMO Synthesis');

    [HMO_Network,RxnCheck] = addReaction(HMO_Network,'HMOSecr_FLNH',['FLNH;HMO[c] ->'],[],false,0,1000,0,'HMO Secretion');

    [HMO_Network,RxnCheck] = addReaction(HMO_Network,'DFLNH1',{'Fa2Ab3GNb3(Ab4(Fa3)GNb6)Ab4G;HMO[c]','DFLNH;HMO[c]'},[-1 1],false,0,1000,0,'HMO Synthesis');
    [HMO_Network,RxnCheck] = addReaction(HMO_Network,'DFLNH2',{'Ab3(Fa4)GNb3(Ab4(Fa3)GNb6)Ab4G;HMO[c]','DFLNH;HMO[c]'},[-1 1],false,0,1000,0,'HMO Synthesis');
    [HMO_Network,RxnCheck] = addReaction(HMO_Network,'DFLNH3',{'Ab3(Fa4)GNb3Ab4(Fa3)GNb3Ab4G;HMO[c]','DFLNH;HMO[c]'},[-1 1],false,0,1000,0,'HMO Synthesis');
    [HMO_Network,RxnCheck] = addReaction(HMO_Network,'DFLNH4',{'Fa2Ab3GNb3Ab4(Fa3)GNb3Ab4G;HMO[c]','DFLNH;HMO[c]'},[-1 1],false,0,1000,0,'HMO Synthesis');
    [HMO_Network,RxnCheck] = addReaction(HMO_Network,'DFLNH5',{'Ab4(Fa3)GNb3Ab4(Fa3)GNb3Ab4G;HMO[c]','DFLNH;HMO[c]'},[-1 1],false,0,1000,0,'HMO Synthesis');
    [HMO_Network,RxnCheck] = addReaction(HMO_Network,'DFLNH6',{'Ab4(Fa3)GNb3(Ab4(Fa3)GNb6)Ab4G;HMO[c]','DFLNH;HMO[c]'},[-1 1],false,0,1000,0,'HMO Synthesis');
    [HMO_Network,RxnCheck] = addReaction(HMO_Network,'DFLNH7',{'Fa2Ab3(Fa4)GNb3(Ab4GNb6)Ab4G;HMO[c]','DFLNH;HMO[c]'},[-1 1],false,0,1000,0,'HMO Synthesis');
    
    
    [HMO_Network,RxnCheck] = addReaction(HMO_Network,'HMOSecr_DFLNH',['DFLNH;HMO[c] ->'],[],false,0,1000,0,'HMO Secretion');

    [HMO_Network,RxnCheck] = addReaction(HMO_Network,'FDSLNH1',{'NNa3Ab3(NNa6)GNb3(Fa2Ab4GNb6)Ab4G;HMO[c]','FDSLNH;HMO[c]'},[-1 1],false,0,1000,0,'HMO Synthesis');
    [HMO_Network,RxnCheck] = addReaction(HMO_Network,'FDSLNH2',{'NNa3Ab3(NNa6)GNb3(Ab4(Fa3)GNb6)Ab4G;HMO[c]','FDSLNH;HMO[c]'},[-1 1],false,0,1000,0,'HMO Synthesis');
    [HMO_Network,RxnCheck] = addReaction(HMO_Network,'FDSLNH3',{'NNa3Ab3(Fa4)GNb3(NNa6Ab4GNb6)Ab4G;HMO[c]','FDSLNH;HMO[c]'},[-1 1],false,0,1000,0,'HMO Synthesis');
    [HMO_Network,RxnCheck] = addReaction(HMO_Network,'FDSLNH4',{'NNa6Ab4(Fa3)GNb3(NNa6Ab4GNb6)Ab4G;HMO[c]','FDSLNH;HMO[c]'},[-1 1],false,0,1000,0,'HMO Synthesis');
    [HMO_Network,RxnCheck] = addReaction(HMO_Network,'FDSLNH5',{'NNa6Ab4(Fa3)GNb3(NNa3Ab4GNb6)Ab4G;HMO[c]','FDSLNH;HMO[c]'},[-1 1],false,0,1000,0,'HMO Synthesis');
    [HMO_Network,RxnCheck] = addReaction(HMO_Network,'FDSLNH6',{'NNa3Ab4(Fa3)GNb3(NNa6Ab4GNb6)Ab4G;HMO[c]','FDSLNH;HMO[c]'},[-1 1],false,0,1000,0,'HMO Synthesis');
    [HMO_Network,RxnCheck] = addReaction(HMO_Network,'FDSLNH7',{'NNa3Ab4(Fa3)GNb3(NNa3Ab4GNb6)Ab4G;HMO[c]','FDSLNH;HMO[c]'},[-1 1],false,0,1000,0,'HMO Synthesis');

    [HMO_Network,RxnCheck] = addReaction(HMO_Network,'HMOSecr_FDSLNH',['FDSLNH;HMO[c] ->'],[],false,0,1000,0,'HMO Secretion');

    [HMO_Network,RxnCheck] = addReaction(HMO_Network,'DSLNH1',{'NNa3Ab3GNb3(NNa6Ab4GNb6)Ab4G;HMO[c]','DSLNH;HMO[c]'},[-1 1],false,0,1000,0,'HMO Synthesis');
    [HMO_Network,RxnCheck] = addReaction(HMO_Network,'DSLNH2',{'NNa3Ab3(NNa6)GNb3(Ab4GNb6)Ab4G;HMO[c]','DSLNH;HMO[c]'},[-1 1],false,0,1000,0,'HMO Synthesis');
    [HMO_Network,RxnCheck] = addReaction(HMO_Network,'DSLNH3',{'NNa6Ab4GNb3(NNa6Ab4GNb6)Ab4G;HMO[c]','DSLNH;HMO[c]'},[-1 1],false,0,1000,0,'HMO Synthesis');

    [HMO_Network,RxnCheck] = addReaction(HMO_Network,'HMOSecr_DSLNH',['DSLNH;HMO[c] ->'],[],false,0,1000,0,'HMO Secretion');
   pause
    
%% STEP5 - Generation of a specific network for the 16 measured HMOs
    nCPU_fva = 4;   
    nSamples = 1e4;
    nThreads = 4;
    tol = 1e-6;
    nSteps=1e5;
    % Set all HMO secretions to zero except for the ones in WT_Profile
    FullNetwork=HMO_Network;
    HMOSecretionIDs = find(strcmp(FullNetwork.subSystems,'HMO Secretion'));

    FullNetwork.ub(HMOSecretionIDs) = 0;
    AuxList = cell(length(WT_Profile.List),1);
    for i=1:length(WT_Profile.List)
        AuxList{i} = [WT_Profile.List{i} '[c]'];
    end;
    ProfileIDs = findSinkIDs(AuxList,FullNetwork,'HMO Secretion');
    FullNetwork.lb(ProfileIDs) = WT_Profile.Fractions'; 
    FullNetwork.ub(ProfileIDs) = WT_Profile.Fractions';
    FullNetwork.description = 'Full HMO Network';

    Network_Red = fastReduceModel(FullNetwork,tol,nCPU_fva,false,false);
    Network_Red.description = 'Reduced HMO Network';
