function [NewNetwork,RxnN] = NetworkEnlarge(Network,RxnFormula,NewGlycan,RxnN,enzyme)

% Addition of synthesis and secretion reactions for HMO network generation


RxnN.syn = RxnN.syn + 1;
NewRxnName = ['HMOSyn_' num2str(RxnN.syn)];
[NewNetwork,RxnCheck] = addReaction(Network,NewRxnName,RxnFormula,[],false,0,1000,0,'HMO Synthesis');
if ~isempty(RxnCheck)
    disp('*** ERROR; DUPLICATE SYNTHESIS REACTIONS !!! ***');
    pause;
end;
if ~ismember(enzyme,NewNetwork.genes)
    NewNetwork.genes{end+1} = enzyme;
    NewNetwork.geneLoc{end+1} = 'c';
end;
geneID = findGeneIDs(NewNetwork,enzyme);
NewNetwork.rxnGeneMat(end,geneID) = true;         % /// addReaction adds a new zero line in rxnGeneMat


if ~hasDMRxn(Network,NewGlycan)
    RxnN.sec = RxnN.sec + 1;
    NewSecrRxnName = ['HMOSecr_' num2str(RxnN.sec)];
    NewSecrRxnFormula = [NewGlycan ' -> '];
    [NewNetwork,RxnCheck] = addReaction(NewNetwork,NewSecrRxnName,NewSecrRxnFormula,[],false,0,1000,0,'HMO Secretion');
    if ~isempty(RxnCheck)
        disp('*** ERROR; DUPLICATE SECRETION REACTIONS !!! ***');
        pause;
    end;                
    NewNetwork.rxnGeneMat(end,1) = true;
end;
            

        
end
