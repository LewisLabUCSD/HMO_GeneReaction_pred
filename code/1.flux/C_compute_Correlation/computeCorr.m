function [P1_data1 P2_data1 P3_data1 P4_data1 P5_data1 P6_data1 P7_data1 P8_data1 P9_data1 P10_data1 ...
          P1_data1_secretor P2_data1_secretor P3_data1_secretor P4_data1_secretor P5_data1_secretor P6_data1_secretor P7_data1_secretor P8_data1_secretor P9_data1_secretor P10_data1_secretor... 
          P1_data2 P2_data2 P3_data2 P4_data2 P5_data2 P6_data2 P7_data2 P8_data2 P9_data2 P10_data2... 
          P1_data2_secretor P2_data2_secretor P3_data2_secretor P4_data2_secretor P5_data2_secretor P6_data2_secretor P7_data2_secretor P8_data2_secretor P9_data2_secretor P10_data2_secretor ...
          R1_data1 R2_data1 R3_data1 R4_data1 R5_data1 R6_data1 R7_data1 R8_data1 R9_data1 R10_data1 ...
          R1_data1_secretor R2_data1_secretor R3_data1_secretor R4_data1_secretor R5_data1_secretor R6_data1_secretor R7_data1_secretor R8_data1_secretor R9_data1_secretor R10_data1_secretor ...
          R1_data2 R2_data2 R3_data2 R4_data2 R5_data2 R6_data2 R7_data2 R8_data2 R9_data2 R10_data2 ...
          R1_data2_secretor R2_data2_secretor R3_data2_secretor R4_data2_secretor R5_data2_secretor R6_data2_secretor R7_data2_secretor R8_data2_secretor R9_data2_secretor R10_data2_secretor model_impossible] = computeCorr(ID_model,solution,Network_Red,normalize_flux)

P1_data1=[];P2_data1=[];P3_data1=[];P4_data1=[];P5_data1=[];P6_data1=[];P7_data1=[];P8_data1=[];P9_data1=[];P10_data1=[];
R1_data1=[];R2_data1=[];R3_data1=[];R4_data1=[];R5_data1=[];R6_data1=[];R7_data1=[];R8_data1=[];R9_data1=[];R10_data1=[];
P1_data1_secretor=[];P2_data1_secretor=[];P3_data1_secretor=[];P4_data1_secretor=[];P5_data1_secretor=[];P6_data1_secretor=[];P7_data1_secretor=[];P8_data1_secretor=[];P9_data1_secretor=[];P10_data1_secretor=[];
R1_data1_secretor=[];R2_data1_secretor=[];R3_data1_secretor=[];R4_data1_secretor=[];R5_data1_secretor=[];R6_data1_secretor=[];R7_data1_secretor=[];R8_data1_secretor=[];R9_data1_secretor=[];R10_data1_secretor=[];
P1_data2=[];P2_data2=[];P3_data2=[];P4_data2=[];P5_data2=[];P6_data2=[];P7_data2=[];P8_data2=[];P9_data2=[];P10_data2=[];
R1_data2=[];R2_data2=[];R3_data2=[];R4_data2=[];R5_data2=[];R6_data2=[];R7_data2=[];R8_data2=[];R9_data2=[];R10_data2=[];
P1_data2_secretor=[];P2_data2_secretor=[];P3_data2_secretor=[];P4_data2_secretor=[];P5_data2_secretor=[];P6_data2_secretor=[];P7_data2_secretor=[];P8_data2_secretor=[];P9_data2_secretor=[];P10_data2_secretor=[];
R1_data2_secretor=[];R2_data2_secretor=[];R3_data2_secretor=[];R4_data2_secretor=[];R5_data2_secretor=[];R6_data2_secretor=[];R7_data2_secretor=[];R8_data2_secretor=[];R9_data2_secretor=[];R10_data2_secretor=[];


N_A ={'LALBA'};
b3GnT_gene ={'B3GNT2','B3GNT3','B3GNT4','B3GNT8','B3GNTL1'}; %L1	 	
a2FucT_gene	={'POFUT1b','POFUT2a','POFUT2b','FUT11','FUT2','FUT3','FUT4','FUT6b','FUT8'};%L2
a3FucT_gene	={'POFUT1b','POFUT2a','POFUT2b','FUT11','FUT2','FUT3','FUT4','FUT6b','FUT8'};%L3
ST3GalT_gene ={'ST3GAL1a','ST3GAL2','ST3GAL3b','ST3GAL4','ST3GAL5b'};%L4				
ST6Gal_gene	={'ST6GAL1b'};						%L5	
b3GalT_gene	={'B3GALNT2','B3GALT4','B3GALT5','B3GALT6'};	%L6	
b4GalT_gene	={'B4GALNT4','B4GALT1','B4GALT2a','B4GALT3','B4GALT4a','B4GALT6','B4GALT7'};%L7
b6GnT_gene ={'GCNT1','GCNT2c','GCNT3'};	 		%L8	
a4FucT_gene	={'POFUT1b','POFUT2a','POFUT2b','FUT11','FUT2','FUT3','FUT4','FUT6b','FUT8'};%L9
ST6GnT_gene	={'ST6GALNAC1','ST6GALNAC2','ST6GALNAC4a','ST6GALNAC6'};%L10

linkage_gene={N_A b3GnT_gene a2FucT_gene a3FucT_gene ST3GalT_gene  ST6Gal_gene b3GalT_gene	 b4GalT_gene	 b6GnT_gene  a4FucT_gene ST6GnT_gene};	 							
load('HMO_data');
%HMO = HMO1;
linkage_geneID={};
for i=1:length(linkage_gene)
	link=linkage_gene{i};
	find_gene=[];
	for j=1:length(link)
        find_gene=[find_gene find(strcmp(link{j},HMO.gene_name)==1)];
    end
	linkage_geneID_data1{i}=find_gene;
end
load('HMO_data2');    
%HMO = HMO2;
linkage_geneID={};
for i=1:length(linkage_gene)
	link=linkage_gene{i};
	find_gene=[];
	for j=1:length(link)
        find_gene=[find_gene find(strcmp(link{j},HMO.gene_name)==1)];
    end
	linkage_geneID_data2{i}=find_gene;
end

model_impossible=[];


%% Loop that compute the correlation for the different dataset for each
%% model defined by ID_model


    display (['Model evaluated = ',num2str(ID_model)])
    %% Generate the minimal network model from solution of MILP in a List of
    %% model structure

    has0=find(solution==0);
    model=removeRxns(Network_Red,Network_Red.rxns(has0)); 

    %% Check that the model is associated with a unique flux distribution
    List_present={'HMO_Init','HMOSyn_2','HMOSecr_2','HMOSecr_3','HMOSyn_4','HMOSecr_4','HMOSecr_6','HMOSecr_7','HMOSecr_15',...
            'HMOSecr_17','HMOSecr_20','HMOSecr_24','HMOSecr_28','HMOSecr_58','HMOSecr_DFLNT','HMOSecr_FLNH','HMOSecr_DFLNH','HMOSecr_FDSLNH','HMOSecr_DSLNH'};
    Present=findRxnIDs(model,List_present);

    Current=1:length(model.rxns);
    RxnEval= Current(~ismember(Current, Present));
    maxFlux=[];
    minFlux=[];

    options = optimset('Display','off','TolFun',1e-6);
    for k=1:length(RxnEval)
        model.c=zeros(length(model.rxns),1);
        model.c(RxnEval(k))=-1;
%         [x,fval] = linprog(model.c,full(model.S),model.b,[],[],model.lb,model.ub,[],options);
%         maxFlux(k)=-fval;
%         model.c(RxnEval(k))=1;
%         [x,fval] = linprog(model.c,full(model.S),model.b,[],[],model.lb,model.ub,[],options);
%         minFlux(k)=fval;
        
        FBAsol = optimizeCbModel(model,'max');
        maxFlux(k)=FBAsol.f;
        FBAsol = optimizeCbModel(model,'min');
        minFlux(k)=FBAsol.f;
    end 
    minFlux(minFlux<1e-6)=0;
    maxFlux(maxFlux<1e-6)=0;
    
    if sum(roundn(minFlux,-4)~=roundn(maxFlux,-4))~=0
        model_impossible= ID_model;
            P1_data1=ones(length(linkage_geneID_data1{2}),1);
            P2_data1=ones(length(linkage_geneID_data1{3}),1);
            P3_data1=ones(length(linkage_geneID_data1{4}),1);
            P4_data1=ones(length(linkage_geneID_data1{5}),1);
            P5_data1=ones(length(linkage_geneID_data1{6}),1);
            P6_data1=ones(length(linkage_geneID_data1{7}),1);
            P7_data1=ones(length(linkage_geneID_data1{8}),1);
            P8_data1=ones(length(linkage_geneID_data1{9}),1);
            P9_data1=ones(length(linkage_geneID_data1{10}),1);
            P10_data1=ones(length(linkage_geneID_data1{11}),1);

            R1_data1=zeros(length(linkage_geneID_data1{2}),1);
            R2_data1=zeros(length(linkage_geneID_data1{3}),1);
            R3_data1=zeros(length(linkage_geneID_data1{4}),1);
            R4_data1=zeros(length(linkage_geneID_data1{5}),1);
            R5_data1=zeros(length(linkage_geneID_data1{6}),1);
            R6_data1=zeros(length(linkage_geneID_data1{7}),1);
            R7_data1=zeros(length(linkage_geneID_data1{8}),1);
            R8_data1=zeros(length(linkage_geneID_data1{9}),1);
            R9_data1=zeros(length(linkage_geneID_data1{10}),1);
            R10_data1=zeros(length(linkage_geneID_data1{11}),1);

            P1_data1_secretor=ones(length(linkage_geneID_data1{2}),1);
            P2_data1_secretor=ones(length(linkage_geneID_data1{3}),1);
            P3_data1_secretor=ones(length(linkage_geneID_data1{4}),1);
            P4_data1_secretor=ones(length(linkage_geneID_data1{5}),1);
            P5_data1_secretor=ones(length(linkage_geneID_data1{6}),1);
            P6_data1_secretor=ones(length(linkage_geneID_data1{7}),1);
            P7_data1_secretor=ones(length(linkage_geneID_data1{8}),1);
            P8_data1_secretor=ones(length(linkage_geneID_data1{9}),1);
            P9_data1_secretor=ones(length(linkage_geneID_data1{10}),1);
            P10_data1_secretor=ones(length(linkage_geneID_data1{11}),1);

            R1_data1_secretor=zeros(length(linkage_geneID_data1{2}),1);
            R2_data1_secretor=zeros(length(linkage_geneID_data1{3}),1);
            R3_data1_secretor=zeros(length(linkage_geneID_data1{4}),1);
            R4_data1_secretor=zeros(length(linkage_geneID_data1{5}),1);
            R5_data1_secretor=zeros(length(linkage_geneID_data1{6}),1);
            R6_data1_secretor=zeros(length(linkage_geneID_data1{7}),1);
            R7_data1_secretor=zeros(length(linkage_geneID_data1{8}),1);
            R8_data1_secretor=zeros(length(linkage_geneID_data1{9}),1);
            R9_data1_secretor=zeros(length(linkage_geneID_data1{10}),1);
            R10_data1_secretor=zeros(length(linkage_geneID_data1{11}),1);
    
            P1_data2=ones(length(linkage_geneID_data2{2}),1);
            P2_data2=ones(length(linkage_geneID_data2{3}),1);
            P3_data2=ones(length(linkage_geneID_data2{4}),1);
            P4_data2=ones(length(linkage_geneID_data2{5}),1);
            P5_data2=ones(length(linkage_geneID_data2{6}),1);
            P6_data2=ones(length(linkage_geneID_data2{7}),1);
            P7_data2=ones(length(linkage_geneID_data2{8}),1);
            P8_data2=ones(length(linkage_geneID_data2{9}),1);
            P9_data2=ones(length(linkage_geneID_data2{10}),1);
            P10_data2=ones(length(linkage_geneID_data2{11}),1);

            R1_data2=zeros(length(linkage_geneID_data2{2}),1);
            R2_data2=zeros(length(linkage_geneID_data2{3}),1);
            R3_data2=zeros(length(linkage_geneID_data2{4}),1);
            R4_data2=zeros(length(linkage_geneID_data2{5}),1);
            R5_data2=zeros(length(linkage_geneID_data2{6}),1);
            R6_data2=zeros(length(linkage_geneID_data2{7}),1);
            R7_data2=zeros(length(linkage_geneID_data2{8}),1);
            R8_data2=zeros(length(linkage_geneID_data2{9}),1);
            R9_data2=zeros(length(linkage_geneID_data2{10}),1);
            R10_data2=zeros(length(linkage_geneID_data2{11}),1);


            P1_data2_secretor=ones(length(linkage_geneID_data2{2}),1);
            P2_data2_secretor=ones(length(linkage_geneID_data2{3}),1);
            P3_data2_secretor=ones(length(linkage_geneID_data2{4}),1);
            P4_data2_secretor=ones(length(linkage_geneID_data2{5}),1);
            P5_data2_secretor=ones(length(linkage_geneID_data2{6}),1);
            P6_data2_secretor=ones(length(linkage_geneID_data2{7}),1);
            P7_data2_secretor=ones(length(linkage_geneID_data2{8}),1);
            P8_data2_secretor=ones(length(linkage_geneID_data2{9}),1);
            P9_data2_secretor=ones(length(linkage_geneID_data2{10}),1);
            P10_data2_secretor=ones(length(linkage_geneID_data2{11}),1);

            R1_data2_secretor=zeros(length(linkage_geneID_data2{2}),1);
            R2_data2_secretor=zeros(length(linkage_geneID_data2{3}),1);
            R3_data2_secretor=zeros(length(linkage_geneID_data2{4}),1);
            R4_data2_secretor=zeros(length(linkage_geneID_data2{5}),1);
            R5_data2_secretor=zeros(length(linkage_geneID_data2{6}),1);
            R6_data2_secretor=zeros(length(linkage_geneID_data2{7}),1);
            R7_data2_secretor=zeros(length(linkage_geneID_data2{8}),1);
            R8_data2_secretor=zeros(length(linkage_geneID_data2{9}),1);
            R9_data2_secretor=zeros(length(linkage_geneID_data2{10}),1);
            R10_data2_secretor=zeros(length(linkage_geneID_data2{11}),1);
    else


    %%  set up constraints
    model.ub=10*ones(length(model.rxns),1);
    model.lb=zeros(length(model.rxns),1);
    model.lb(1)=1;
    model.ub(1)=1;

    %% Compute the coefficient related to each reaction for further sum in
    %% correlation construction
        A= full(model.S);
        B = sum(A,1);
        coeff_flux={};
        for i=1:length(B)
            if B(i)==1
                coeff_flux{i}= 1;
            elseif B(i)==-1
                coeff_flux{i}= i;
            elseif B(i)==0
                nb_rxn=find(A(:,i)==-1);
                nb=find(A(nb_rxn,:)==1);
                coeff_flux{i}= nb;
            end
        end

    %% Find the ID of HMO secretion reaction
    % 2FL HMOSecr_2	Fa2Ab4G;HMO [c] 	->	
    % 3FL HMOSecr_3	Ab4(Fa3)G;HMO[c] 	->	
    % LNnT HMOSecr_7	Ab4GNb3Ab4G;HMO[c] 	->	
    % 3SL HMOSecr_4	NNa3Ab4G;HMO[c] 	->	
    % LNT HMOSecr_6	Ab3GNb3Ab4G;HMO[c] 	->	
    % LNFPI HMOSecr_15	Fa2Ab3GNb3Ab4G;HMO[c] 	->	
    % LNFPII HMOSecr_17	Ab3(Fa4)GNb3Ab4G;HMO[c] 	->	
    % LNFPIII HMOSecr_24	Ab4(Fa3)GNb3Ab4G;HMO[c] 	->	
    % LSTb  HMOSecr_20	Ab3(NNa6)GNb3Ab4G;HMO[c] 	->	
    % LSTc HMOSecr_28	NNa6Ab4GNb3Ab4G;HMO[c] 	->	
    % DSLNT HMOSecr_58	NNa3Ab3(NNa6)GNb3Ab4G;HMO[c] 	->	

    List_secr={'HMOSecr_2','HMOSecr_3','HMOSecr_7','HMOSecr_4','HMOSecr_6','HMOSecr_15','HMOSecr_17','HMOSecr_24','HMOSecr_20',...
        'HMOSecr_28','HMOSecr_DFLNT','HMOSecr_58','HMOSecr_FLNH','HMOSecr_DFLNH','HMOSecr_FDSLNH','HMOSecr_DSLNH'};
    Rxn_secr =findRxnIDs(model,List_secr);

    %% Compute FBA solution for the model for data 1
    load('HMO_data');               
    %HMO = HMO1;
    O_exp_rel=[];
    for j =1:length(HMO.conc(1,:)) %47
        for i =1:length(HMO.conc(:,1)) %17
        O_exp_rel(i,j)=(HMO.conc(i,j)./sum(HMO.conc(:,j))); 
        end
    end
    HMO.conc_Fractions=[O_exp_rel];

    FluxSol=[];

    for j = 1:length(HMO.conc_Fractions(1,:))
        model.lb(Rxn_secr) = HMO.conc_Fractions(:,j); 
        model.ub(Rxn_secr) = HMO.conc_Fractions(:,j);
%         mins=linprog(model.c,full(model.S),model.b,[],[],model.lb,model.ub,[],options);
%         mins(mins<1e-6)=0;
%         FluxSol=[FluxSol mins];
        
        FBAsol = optimizeCbModel(model,'min');
        mins=FBAsol.x;
        FluxSol=[FluxSol mins];
    end


    %% Compute correlation for all potential genes-linkage combinations by
    %% using an extrapolation of flux values for all mRNA sample time and using
    %% the normalized flux values

    %% DATA1 - ALL
    G_name=HMO.gene_name;
    G_exp= HMO.gene(:,find(ismember(HMO.gene_label,HMO.conc_label)==1));

    R_model={};
    P_model={};
    FluxRatio={};

    flux=FluxSol'; %need to check the dimension '   

    if normalize_flux
        flux_norm=[];
        for i = 1:length(flux(1,:))
            if coeff_flux{i}==1
                if sum(flux(:,i))==0
                    flux_norm(:,i) =flux(:,i);
                else
                    flux_norm(:,i) =flux(:,i)./flux(:,coeff_flux{i});
                end
            else
                for j=1: length(coeff_flux{i})
                    if sum(flux(:,coeff_flux{i}(j)))==0
                        j=j+1;
                    else
                        flux_norm(:,i) =flux(:,i)./flux(:,coeff_flux{i}(j));
                    end
                end
            end
        end
    else
        flux_norm=flux;
    end

    R_score=zeros(9,10);
    P_score=ones(9,10);
    [row,col]=find(model.rxnGeneMat(:,2:end)==1);
    %due to model.rxnGeneMat(:,2:end) - we not taken into account N/A
    %related rxns
    col=col+1;
    score.rxnsID=model.rxns(row);
    score.linkage=model.genes(col);
    
    flux_sum=[];
    for n=2:11
        flux_sum=[flux_sum sum(flux_norm(:,row(find(col==n))),2)];
    end

    for j= 1:10

        for i = 1:length(linkage_geneID_data1{j+1})

            G_ID=linkage_geneID_data1{j+1};
            G_mes= G_exp(G_ID(i),:);

            [corr_value,p_val] = corr(flux_sum(:,j),(G_mes)','type','Spearman');
            R_score(i,j) = (corr_value);
            P_score(i,j) = (p_val);
        end       
    end

    R_score(isnan(R_score))=0;
    P_score(isnan(P_score))=1;
    FluxRatio=flux_norm;
    

    P1_data1=P_score(1:length(linkage_geneID_data1{2}),1);
    P2_data1=P_score(1:length(linkage_geneID_data1{3}),2);
    P3_data1=P_score(1:length(linkage_geneID_data1{4}),3);
    P4_data1=P_score(1:length(linkage_geneID_data1{5}),4);
    P5_data1=P_score(1:length(linkage_geneID_data1{6}),5);
    P6_data1=P_score(1:length(linkage_geneID_data1{7}),6);
    P7_data1=P_score(1:length(linkage_geneID_data1{8}),7);
    P8_data1=P_score(1:length(linkage_geneID_data1{9}),8);
    P9_data1=P_score(1:length(linkage_geneID_data1{10}),9);
    P10_data1=P_score(1:length(linkage_geneID_data1{11}),10);

    R1_data1=R_score(1:length(linkage_geneID_data1{2}),1);
    R2_data1=R_score(1:length(linkage_geneID_data1{3}),2);
    R3_data1=R_score(1:length(linkage_geneID_data1{4}),3);
    R4_data1=R_score(1:length(linkage_geneID_data1{5}),4);
    R5_data1=R_score(1:length(linkage_geneID_data1{6}),5);
    R6_data1=R_score(1:length(linkage_geneID_data1{7}),6);
    R7_data1=R_score(1:length(linkage_geneID_data1{8}),7);
    R8_data1=R_score(1:length(linkage_geneID_data1{9}),8);
    R9_data1=R_score(1:length(linkage_geneID_data1{10}),9);
    R10_data1=R_score(1:length(linkage_geneID_data1{11}),10);

    %% DATA1 - SECRETOR
    G_name=HMO.gene_name;
    G_exp= HMO.gene(:,find(ismember(HMO.gene_label,HMO.conc_label)==1));
    G_exp=[G_exp(:,1:8) G_exp(:,25:end)];

    R_model={};
    P_model={};
    FluxRatio={};
    
    flux=FluxSol';
    flux=[flux(1:8,:); flux(25:end,:)];    

    flux_norm=[];
    for i = 1:length(flux(1,:))

        if coeff_flux{i}==1
            if sum(flux(:,i))==0
                flux_norm(:,i) =flux(:,i);
            else
                flux_norm(:,i) =flux(:,i)./flux(:,coeff_flux{i});
            end
        else
            for j=1: length(coeff_flux{i})
                if sum(flux(:,coeff_flux{i}(j)))==0
                    j=j+1;
                else
                    flux_norm(:,i) =flux(:,i)./flux(:,coeff_flux{i}(j));
                end
            end
        end
    end

    R_score=zeros(9,10);
    P_score=ones(9,10);
    [row,col]=find(model.rxnGeneMat(:,2:end)==1);
    %due to model.rxnGeneMat(:,2:end) - we not taken into account N/A
    %related rxns
    col=col+1;
    score.rxnsID=model.rxns(row);
    score.linkage=model.genes(col);
    flux_sum=[];
    for n=2:11
        flux_sum=[flux_sum sum(flux_norm(:,row(find(col==n))),2)];
    end
    
    for j= 1:10
        for i = 1:length(linkage_geneID_data1{j+1})
            G_ID=linkage_geneID_data1{j+1};
            G_mes= G_exp(G_ID(i),:);
            [corr_value,p_val] = corr(flux_sum(:,j),(G_mes)','type','Spearman');
            R_score(i,j) = (corr_value);
            P_score(i,j) = (p_val);
        end       
    end
    R_score(isnan(R_score))=0;
    P_score(isnan(P_score))=1;

    
    P1_data1_secretor=P_score(1:length(linkage_geneID_data1{2}),1);
    P2_data1_secretor=P_score(1:length(linkage_geneID_data1{3}),2);
    P3_data1_secretor=P_score(1:length(linkage_geneID_data1{4}),3);
    P4_data1_secretor=P_score(1:length(linkage_geneID_data1{5}),4);
    P5_data1_secretor=P_score(1:length(linkage_geneID_data1{6}),5);
    P6_data1_secretor=P_score(1:length(linkage_geneID_data1{7}),6);
    P7_data1_secretor=P_score(1:length(linkage_geneID_data1{8}),7);
    P8_data1_secretor=P_score(1:length(linkage_geneID_data1{9}),8);
    P9_data1_secretor=P_score(1:length(linkage_geneID_data1{10}),9);
    P10_data1_secretor=P_score(1:length(linkage_geneID_data1{11}),10);

    R1_data1_secretor=R_score(1:length(linkage_geneID_data1{2}),1);
    R2_data1_secretor=R_score(1:length(linkage_geneID_data1{3}),2);
    R3_data1_secretor=R_score(1:length(linkage_geneID_data1{4}),3);
    R4_data1_secretor=R_score(1:length(linkage_geneID_data1{5}),4);
    R5_data1_secretor=R_score(1:length(linkage_geneID_data1{6}),5);
    R6_data1_secretor=R_score(1:length(linkage_geneID_data1{7}),6);
    R7_data1_secretor=R_score(1:length(linkage_geneID_data1{8}),7);
    R8_data1_secretor=R_score(1:length(linkage_geneID_data1{9}),8);
    R9_data1_secretor=R_score(1:length(linkage_geneID_data1{10}),9);
    R10_data1_secretor=R_score(1:length(linkage_geneID_data1{11}),10);


    %% Compute FBA solution for the model for data 2
    load('HMO_data2');               
    %HMO = HMO2;
    O_exp_rel=[];
    for j =1:length(HMO.conc(1,:)) %10
        for i =1:length(HMO.conc(:,1)) %17
        O_exp_rel(i,j)=(HMO.conc(i,j)./sum(HMO.conc(:,j))); 
        end
    end
    HMO.conc_Fractions=[O_exp_rel];

    FluxSol=[];

    for j = 1:length(HMO.conc_Fractions(1,:))
        model.lb(Rxn_secr) = HMO.conc_Fractions(:,j); 
        model.ub(Rxn_secr) = HMO.conc_Fractions(:,j);
%         mins=linprog(model.c,full(model.S),model.b,[],[],model.lb,model.ub,[],options);
%         mins(mins<1e-6)=0;
%         FluxSol=[FluxSol mins];
        FBAsol = optimizeCbModel(model,'min');
        mins=FBAsol.x;
        FluxSol=[FluxSol mins];
    end


    %% Compute correlation for all potential genes-linkage combinations by
    %% using an extrapolation of flux values for all mRNA sample time and using
    %% the normalized flux values

    %% DATA2 - ALL
    G_name=HMO.gene_name;
    G_exp= HMO.gene(:,find(ismember(HMO.gene_label,HMO.conc_label)==1));

    R_model={};
    P_model={};
    FluxRatio={};

    flux=FluxSol'; %need to check the dimension '  

    if normalize_flux
        flux_norm=[];
        for i = 1:length(flux(1,:))
            if coeff_flux{i}==1
                if sum(flux(:,i))==0
                    flux_norm(:,i) =flux(:,i);
                else
                    flux_norm(:,i) =flux(:,i)./flux(:,coeff_flux{i});
                end
            else
                for j=1: length(coeff_flux{i})
                    if sum(flux(:,coeff_flux{i}(j)))==0
                        j=j+1;
                    else
                        flux_norm(:,i) =flux(:,i)./flux(:,coeff_flux{i}(j));
                    end
                end
            end
        end
    else
        flux_norm=flux;
    end

    R_score=zeros(9,10);
    P_score=ones(9,10);
    [row,col]=find(model.rxnGeneMat(:,2:end)==1);
    %due to model.rxnGeneMat(:,2:end) - we not taken into account N/A
    %related rxns
    col=col+1;
    score.rxnsID=model.rxns(row);
    score.linkage=model.genes(col);
    flux_sum=[];
    
    for n=2:11
        flux_sum=[flux_sum sum(flux_norm(:,row(find(col==n))),2)];
    end
    
    for j= 1:10
        for i = 1:length(linkage_geneID_data2{j+1})
            G_ID=linkage_geneID_data2{j+1};
            G_mes= G_exp(G_ID(i),:);
            [corr_value,p_val] = corr(flux_sum(:,j),(G_mes)','type','Spearman');
            R_score(i,j) = (corr_value);
            P_score(i,j) = (p_val);
        end       
    end
    R_score(isnan(R_score))=0;
    P_score(isnan(P_score))=1;
    
    FluxRatio=flux_norm;
    

    P1_data2=P_score(1:length(linkage_geneID_data2{2}),1);
    P2_data2=P_score(1:length(linkage_geneID_data2{3}),2);
    P3_data2=P_score(1:length(linkage_geneID_data2{4}),3);
    P4_data2=P_score(1:length(linkage_geneID_data2{5}),4);
    P5_data2=P_score(1:length(linkage_geneID_data2{6}),5);
    P6_data2=P_score(1:length(linkage_geneID_data2{7}),6);
    P7_data2=P_score(1:length(linkage_geneID_data2{8}),7);
    P8_data2=P_score(1:length(linkage_geneID_data2{9}),8);
    P9_data2=P_score(1:length(linkage_geneID_data2{10}),9);
    P10_data2=P_score(1:length(linkage_geneID_data2{11}),10);

    R1_data2=R_score(1:length(linkage_geneID_data2{2}),1);
    R2_data2=R_score(1:length(linkage_geneID_data2{3}),2);
    R3_data2=R_score(1:length(linkage_geneID_data2{4}),3);
    R4_data2=R_score(1:length(linkage_geneID_data2{5}),4);
    R5_data2=R_score(1:length(linkage_geneID_data2{6}),5);
    R6_data2=R_score(1:length(linkage_geneID_data2{7}),6);
    R7_data2=R_score(1:length(linkage_geneID_data2{8}),7);
    R8_data2=R_score(1:length(linkage_geneID_data2{9}),8);
    R9_data2=R_score(1:length(linkage_geneID_data2{10}),9);
    R10_data2=R_score(1:length(linkage_geneID_data2{11}),10);

    %% DATA2 - SECRETOR
    G_name=HMO.gene_name;
    G_exp= HMO.gene(:,find(ismember(HMO.gene_label,HMO.conc_label)==1));
    G_exp=G_exp(:,3:end);

    R_model={};
    P_model={};
    FluxRatio={};
    
    flux=FluxSol';
    flux=[flux(3:end,:)];    

    flux_norm=[];
    for i = 1:length(flux(1,:))

        if coeff_flux{i}==1
            if sum(flux(:,i))==0
                flux_norm(:,i) =flux(:,i);
            else
                flux_norm(:,i) =flux(:,i)./flux(:,coeff_flux{i});
            end
        else
            for j=1: length(coeff_flux{i})
                if sum(flux(:,coeff_flux{i}(j)))==0
                    j=j+1;
                else
                    flux_norm(:,i) =flux(:,i)./flux(:,coeff_flux{i}(j));
                end
            end
        end
    end

    R_score=zeros(9,10);
    P_score=ones(9,10);
    [row,col]=find(model.rxnGeneMat(:,2:end)==1);
    %due to model.rxnGeneMat(:,2:end) - we not taken into account N/A
    %related rxns
    col=col+1;
    score.rxnsID=model.rxns(row);
    score.linkage=model.genes(col);
    flux_sum=[];
    
    for n=2:11
        flux_sum=[flux_sum sum(flux_norm(:,row(find(col==n))),2)];
    end
    
    for j= 1:10
        for i = 1:length(linkage_geneID_data2{j+1})
            G_ID=linkage_geneID_data2{j+1};
            G_mes= G_exp(G_ID(i),:);
            [corr_value,p_val] = corr(flux_sum(:,j),(G_mes)','type','Spearman');
            R_score(i,j) = (corr_value);
            P_score(i,j) = (p_val);
        end       
    end
    R_score(isnan(R_score))=0;
    P_score(isnan(P_score))=1;

    
    P1_data2_secretor=P_score(1:length(linkage_geneID_data2{2}),1);
    P2_data2_secretor=P_score(1:length(linkage_geneID_data2{3}),2);
    P3_data2_secretor=P_score(1:length(linkage_geneID_data2{4}),3);
    P4_data2_secretor=P_score(1:length(linkage_geneID_data2{5}),4);
    P5_data2_secretor=P_score(1:length(linkage_geneID_data2{6}),5);
    P6_data2_secretor=P_score(1:length(linkage_geneID_data2{7}),6);
    P7_data2_secretor=P_score(1:length(linkage_geneID_data2{8}),7);
    P8_data2_secretor=P_score(1:length(linkage_geneID_data2{9}),8);
    P9_data2_secretor=P_score(1:length(linkage_geneID_data2{10}),9);
    P10_data2_secretor=P_score(1:length(linkage_geneID_data2{11}),10);

    R1_data2_secretor=R_score(1:length(linkage_geneID_data2{2}),1);
    R2_data2_secretor=R_score(1:length(linkage_geneID_data2{3}),2);
    R3_data2_secretor=R_score(1:length(linkage_geneID_data2{4}),3);
    R4_data2_secretor=R_score(1:length(linkage_geneID_data2{5}),4);
    R5_data2_secretor=R_score(1:length(linkage_geneID_data2{6}),5);
    R6_data2_secretor=R_score(1:length(linkage_geneID_data2{7}),6);
    R7_data2_secretor=R_score(1:length(linkage_geneID_data2{8}),7);
    R8_data2_secretor=R_score(1:length(linkage_geneID_data2{9}),8);
    R9_data2_secretor=R_score(1:length(linkage_geneID_data2{10}),9);
    R10_data2_secretor=R_score(1:length(linkage_geneID_data2{11}),10);
    end

end

