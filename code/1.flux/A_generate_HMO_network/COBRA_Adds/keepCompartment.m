function newmodel = keepCompartment(model, compartments)
% This function removes reactions in all compartments except those
% specified by the cell array "compartments"
% 
%INPUTS
% model                     COBRA model structure
% compartments              cell array of strings (e.g., to discard all
%                           reactions except those in the mitochondria and 
%                           cytosol, compartments = {'[m]','[c]'};
%
%OUTPUT
% newmodel                  COBRA model with reactions in the specified
%                           compartmetns
% 
% Nathan Lewis / Philipp Spahn
% June 8, 2008 / Aug 6, 2014


% compartments is a cell array list of compartments to keep 
compartments = regexprep(compartments, '\[','\\\[');
compartments = regexprep(compartments, '\]','\\\]');

% make a list of metabolites which are in the desired compartment
Mets2KeepList = {};
k = 0;
for j = 1:length(model.mets)
    RegMatch = regexpi(model.mets{j},compartments);
    if sum(cellfun(@isempty,RegMatch)) < length(compartments)
        k = k+1;
        Mets2KeepList{k} = model.mets{j};
    end;
end;

% make a list of rxns to keep
Rxns2KeepList = findRxnsFromMets(model,Mets2KeepList);

% make a list of rxns to be kicked out
Rxns2KickList = {};
k = 0;
for i = 1:length(model.rxns)
    if ~ismember(model.rxns{i},Rxns2KeepList)
        k = k+1;
        Rxns2KickList{k} = model.rxns{i};
    end;
end;

% kick out all reactions on the Kick-out List
newmodel = removeRxns(model,Rxns2KickList);


end
