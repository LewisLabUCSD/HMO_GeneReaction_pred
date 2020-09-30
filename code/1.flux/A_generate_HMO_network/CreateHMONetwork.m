function HMO_Network = CreateHMONetwork(BasicNetwork,K)

% Generic reaction network of human milk oligosaccharides
%%construction of the generic
% HMO_01 = 'Fa2Ab4G;HMO'; %2FL
% HMO_02 = 'Ab4(Fa3)G;HMO';%3FL
% HMO_03 = 'NNa3Ab4G;HMO'; %3SL
% HMO_04 = 'Ab3GNb3Ab4G;HMO';%LNT
% HMO_05 = 'Ab4GNb3Ab4G;HMO';%LNnT
% HMO_06 = 'Fa2Ab3GNb3Ab4G;HMO';%LNFPI 
% HMO_07 = 'Ab3(Fa4)GNb3Ab4G;HMO';% LNFPII
% HMO_08 = 'Ab4(Fa3)GNb3Ab4G;HMO';% LNFPIII
% HMO_09 = 'Ab3(NNa6)GNb3Ab4G;HMO';% LSTb  
% HMO_10 = 'NNa6Ab4GNb3Ab4G;HMO';% LSTc 
% HMO_11 = 'NNa3Ab3(NNa6)GNb3Ab4G;HMO';%DSLNT 
% HMO_12 = 'Fa2Ab3(Fa4)GNb3Ab4G;HMO';% DFLNT 
% HMO_13 = 'Ab3(Fa4)GNb3Ab4(Fa3)G;HMO';% DFLNT 
% HMO_14 = 'Ab4(Fa3)GNb3Ab4(Fa3)G;HMO';% DFLNT 
% HMO_15 = 'Ab3GNb3Ab4(Fa3)GNb3Ab4G;HMO';% FLNH 	
% HMO_16 = 'Ab3(Fa4)GNb3Ab4GNb3Ab4G;HMO';% FLNH 
% HMO_17 = 'Fa2Ab3GNb3(Ab4GNb6)Ab4G;HMO';% FLNH 	
% HMO_18 = 'Ab3GNb3(Ab4(Fa3)GNb6)Ab4G;HMO';% FLNH 
% HMO_19 = 'Fa2Ab3GNb3(Ab4(Fa3)GNb6)Ab4G;HMO';% DFLNH  
% HMO_20 = 'Ab3(Fa4)GNb3(Ab4(Fa3)GNb6)Ab4G;HMO';% DFLNH  
% HMO_21 = 'Ab3(Fa4)GNb3Ab4(Fa3)GNb3Ab4G;HMO';% DFLNH  	
% HMO_22 = 'Fa2Ab3GNb3Ab4(Fa3)GNb3Ab4G;HMO';% DFLNH  
% HMO_23 = 'Ab4(Fa3)GNb3Ab4(Fa3)GNb3Ab4G;HMO';% DFLNH  	
% HMO_24 = 'Ab4(Fa3)GNb3(Ab4(Fa3)GNb6)Ab4G;HMO';% DFLNH     
% HMO_25 = 'NNa3Ab3(NNa6)GNb3(Fa2Ab4GNb6)Ab4G;HMO';% FDSLNH	
% HMO_26 = 'NNa3Ab3(NNa6)GNb3(Ab4(Fa3)GNb6)Ab4G;HMO';% FDSLNH	
% HMO_27 = 'NNa3Ab3(Fa4)GNb3(NNa6Ab4GNb6)Ab4G;HMO';% FDSLNH	
% HMO_28 = 'NNa6Ab4(Fa3)GNb3(NNa6Ab4GNb6)Ab4G;HMO';% FDSLNH		
% HMO_29 = 'NNa6Ab4(Fa3)GNb3(NNa3Ab4GNb6)Ab4G;HMO';% FDSLNH		
% HMO_30 = 'NNa3Ab4(Fa3)GNb3(NNa6Ab4GNb6)Ab4G;HMO';% FDSLNH	
% HMO_31 = 'NNa3Ab4(Fa3)GNb3(NNa3Ab4GNb6)Ab4G;HMO';% FDSLNH	
% HMO_32 = 'NNa3Ab3GNb3(NNa6Ab4GNb6)Ab4G;HMO';% DSLNH  
% HMO_33 = 'NNa3Ab3(NNa6)GNb3(Ab4GNb6)Ab4G;HMO';% DSLNH  
% HMO_34 = 'NNa6Ab4GNb3(NNa6Ab4GNb6)Ab4G;HMO';% DSLNH  	
%         
% S=G     P=|Ab4G => GAL1
% S=G     P=(Fa3)G => FUC1

% S=Ab     P=|GNb3|Ab => GLCNAC1
% S=Ab     P=GNb6|Ab => GLCNAC2
% S=Ab     P=(Fa2)Ab => FUC2
% S=Ab     P=(NNa3)Ab  => NEUAC1
% S=Ab     P=(NNa6)Ab => NEUAC2

% S=GNb	P=Ab3GNb  => GAL3
% S=GNb	P=Ab4GNb  => GAL4
% S=GNb	P=(Fa3)GNb => FUC3
% S=GNb	P=(Fa4)GNb => FUC4
% S=GNb	P=(NNa6)GNb => NEUAC3


HMO_start = 'Ab4G;HMO[c]';
HMO_List{1} = {HMO_start};
RxnN.syn = 0;
RxnN.sec = 0;

HMO_Network = BasicNetwork;
HMO_Network.description = ['Human Milk Oligosaccharide Network - Complexity Level ' num2str(K)];

k = 1;
H1 = waitbar(0,'Creating Glycan Network ...','Position',[400 400 300 50]);
tic
while k <= K
    NextList = {};
    CurrList = HMO_List{k};
    n = 0;
    H2 = waitbar(0,['Assembling Complexity Level ' num2str(k) ' ... '],'Position',[400 300 300 50]);
    for j = 1:length(CurrList)
        CurrGlycan = CurrList{j};
        %% b3GnT Elongation => GLCNAC1
        %----------------------------------------------------
        SubstrateString = '(A';
        ProductString = '(GNb3A';
        CurrGlyc0 = ['(' CurrGlycan];
        b3GnT = strfind(CurrGlyc0,SubstrateString);
        for i = 1:length(b3GnT)
            f = b3GnT(i);
            NewGlyc0 = [CurrGlyc0(1:(f-1)) ProductString CurrGlyc0((f+2):end)];
            NewGlycan = NewGlyc0(2:end);
            RxnFormula = [CurrGlycan ' -> ' NewGlycan];
            [HMO_Network,RxnN] = NetworkEnlargeHMO(HMO_Network,RxnFormula,NewGlycan,RxnN,'b3GnT');
            if ~ismember(NewGlycan,NextList)
                n = n+1; NextList{n} = NewGlycan;
            end;
        end;
        %% b3GalT Elongation => GAL 3
        %----------------------------------------------------
        SubstrateString = '(GN';
        ProductString = '(Ab3GN';
        CurrGlyc0 = ['(' CurrGlycan];
        b3GalT = strfind(CurrGlyc0,SubstrateString);
        for i = 1:length(b3GalT)
            f = b3GalT(i);
            NewGlyc0 = [CurrGlyc0(1:(f-1)) ProductString CurrGlyc0((f+3):end)];
            NewGlycan = NewGlyc0(2:end);
            RxnFormula = [CurrGlycan ' -> ' NewGlycan];
            [HMO_Network,RxnN] = NetworkEnlargeHMO(HMO_Network,RxnFormula,NewGlycan,RxnN,'b3GalT');
            if ~ismember(NewGlycan,NextList)
                n = n+1; NextList{n} = NewGlycan;
            end;
        end;
        %% b4GalT Elongation => GAL4
        %----------------------------------------------------
        SubstrateString = '(GN';
        ProductString = '(Ab4GN';
        CurrGlyc0 = ['(' CurrGlycan];
        b4GalT = strfind(CurrGlyc0,SubstrateString);
        for i = 1:length(b4GalT)
            f = b4GalT(i);
            NewGlyc0 = [CurrGlyc0(1:(f-1)) ProductString CurrGlyc0((f+3):end)];
            NewGlycan = NewGlyc0(2:end);
            RxnFormula = [ CurrGlycan ' -> ' NewGlycan];
            [HMO_Network,RxnN] = NetworkEnlargeHMO(HMO_Network,RxnFormula,NewGlycan,RxnN,'b4GalT');
            if ~ismember(NewGlycan,NextList)
                n = n+1; NextList{n} = NewGlycan;
            end;
        end;
        %% b6GnT Branching => GLCNAC2
        %---------------------------------------------------- 	
        SubstrateString = 'GNb3Ab4';
        ProductString = 'GNb3(GNb6)Ab4';
        CurrGlyc0 = ['(' CurrGlycan];
        b6GnT = strfind(CurrGlyc0,SubstrateString);
        for i = 1:length(b6GnT)
            f = b6GnT(i);
            NewGlyc0 = [CurrGlyc0(1:(f-1)) ProductString CurrGlyc0((f+7):end)];
            NewGlycan = NewGlyc0(2:end);
            RxnFormula = [CurrGlycan ' -> ' NewGlycan];
            [HMO_Network,RxnN] = NetworkEnlargeHMO(HMO_Network,RxnFormula,NewGlycan,RxnN,'b6GnT');
            if ~ismember(NewGlycan,NextList)
                n = n+1; NextList{n} = NewGlycan;
            end;
        end;
        %% a2FucT Fucosylation => FUC2
        %----------------------------------------------------
        SubstrateString = '(A';
        ProductString = '(Fa2A';
        CurrGlyc0 = ['(' CurrGlycan];
        a2FucT = strfind(CurrGlyc0,SubstrateString);
        for i = 1:length(a2FucT)
            f = a2FucT(i);
            NewGlyc0 = [CurrGlyc0(1:(f-1)) ProductString CurrGlyc0((f+2):end)];
            NewGlycan = NewGlyc0(2:end);
            RxnFormula = [CurrGlycan ' -> ' NewGlycan];
            [HMO_Network,RxnN] = NetworkEnlargeHMO(HMO_Network,RxnFormula,NewGlycan,RxnN,'a2FucT');
            if ~ismember(NewGlycan,NextList)
                n = n+1; NextList{n} = NewGlycan;
            end;
        end;
        %% a3GFucT Fucosylation=> FUC1 & FUC3
        %----------------------------------------------------
        SubstrateString = 'Ab4G';
        ProductString = 'Ab4(Fa3)G';
        CurrGlyc0 = ['(' CurrGlycan];
        a3FucT = strfind(CurrGlyc0,SubstrateString);
        for i = 1:length(a3FucT)
            f = a3FucT(i);
            NewGlyc0 = [CurrGlyc0(1:(f-1)) ProductString CurrGlyc0((f+4):end)];
            NewGlycan = NewGlyc0(2:end);
            RxnFormula = [ CurrGlycan ' -> ' NewGlycan];
            [HMO_Network,RxnN] = NetworkEnlargeHMO(HMO_Network,RxnFormula,NewGlycan,RxnN,'a3FucT');
            if ~ismember(NewGlycan,NextList)
                n = n+1; NextList{n} = NewGlycan;
            end;
        end;

        
        %% a4FucT Fucosylation=> FUC4
        %----------------------------------------------------
        SubstrateString = 'GNb';
        ProductString = '(Fa4)GNb';
        CurrGlyc0 = ['(' CurrGlycan];
        a4FucT = regexp(CurrGlyc0,'(\w[^)]GNb|[^4]\)GNb)');
        for i = 1:length(a4FucT)
            f = a4FucT(i);
            NewGlyc0 = [CurrGlyc0(1:(f+1)) ProductString CurrGlyc0((f+5):end)];
            NewGlycan = NewGlyc0(2:end);
            RxnFormula = [ CurrGlycan ' -> ' NewGlycan];
            [HMO_Network,RxnN] = NetworkEnlargeHMO(HMO_Network,RxnFormula,NewGlycan,RxnN,'a4FucT');
            if ~ismember(NewGlycan,NextList)
                n = n+1; NextList{n} = NewGlycan;
            end;
        end;
        %% ST3GalT Terminal Sialylation => NEUAC1
        %-----------------------------------a6SiaT-----------------
        SubstrateString = '(A';
        ProductString = '(NNa3A';
        CurrGlyc0 = ['(' CurrGlycan];
        a3SiaT = strfind(CurrGlyc0,SubstrateString);
        for i = 1:length(a3SiaT)
            f = a3SiaT(i);
            NewGlyc0 = [CurrGlyc0(1:(f-1)) ProductString CurrGlyc0((f+2):end)];
            NewGlycan = NewGlyc0(2:end);
            RxnFormula = [CurrGlycan ' -> ' NewGlycan];
            [HMO_Network,RxnN] = NetworkEnlargeHMO(HMO_Network,RxnFormula,NewGlycan,RxnN,'ST3Gal');
            if ~ismember(NewGlycan,NextList)
                n = n+1; NextList{n} = NewGlycan;
            end;
        end;    
        %% ST6Gal Terminal Sialylation => NEUAC2
        %----------------------------------------------------
        SubstrateString = '(A';
        ProductString = '(NNa6A';
        CurrGlyc0 = ['(' CurrGlycan];
        a6SiaT = strfind(CurrGlyc0,SubstrateString);
        for i = 1:length(a6SiaT)
            f = a6SiaT(i);
            NewGlyc0 = [CurrGlyc0(1:(f-1)) ProductString CurrGlyc0((f+2):end)];
            NewGlycan = NewGlyc0(2:end);
            RxnFormula = [CurrGlycan ' -> ' NewGlycan];
            [HMO_Network,RxnN] = NetworkEnlargeHMO(HMO_Network,RxnFormula,NewGlycan,RxnN,'ST6Gal');
            if ~ismember(NewGlycan,NextList)
                n = n+1; NextList{n} = NewGlycan;
            end;
        end;         
        %% ST6GnT Subterminal Sialylation => NEUAC3
        %----------------------------------------------------
        SubstrateString = '(Ab3GN';
        ProductString = '(Ab3(NNa6)GN';
        CurrGlyc0 = ['(' CurrGlycan];
        a6SiaT = strfind(CurrGlyc0,SubstrateString);
        for i = 1:length(a6SiaT)
            f = a6SiaT(i);
            NewGlyc0 = [CurrGlyc0(1:(f-1)) ProductString CurrGlyc0((f+6):end)];
            NewGlycan = NewGlyc0(2:end);
            RxnFormula = [ CurrGlycan ' -> ' NewGlycan];
            [HMO_Network,RxnN] = NetworkEnlargeHMO(HMO_Network,RxnFormula,NewGlycan,RxnN,'ST3GnT');
            if ~ismember(NewGlycan,NextList)
                n = n+1; NextList{n} = NewGlycan;
            end;
        end;           
        waitbar(j/length(CurrList),H2);
    end;
    close(H2);
    waitbar(k/K,H1);
    HMO_List{k+1} = NextList';
    k = k+1;   
end;
toc
close(H1);


end