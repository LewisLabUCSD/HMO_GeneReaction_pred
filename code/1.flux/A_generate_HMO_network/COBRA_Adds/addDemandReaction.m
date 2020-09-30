function [model,rxnNames] = addDemandReaction(model,metaboliteNameList,RxnName)
% addDemandReaction adds demand reactions for a set of metabolites
% The reaction names for the demand reactions will be DM_[metaboliteName]
%
% model = addDemandReaction(model,metaboliteNameList)
%
% INPUTS
% model                  COBRA model structure
% metaboliteNameList     List of metabolite names (cell array)
% RxnName (optional)     Reaction name
%
% OUTPUTS
% model                  COBRA model structure with added demand reactions
% rxnNames               List of added reactions
%
% Markus Herrgard 5/8/07
% Ines Thiele 03/09 - Corrected reaction coefficient for demand reaction

if (~iscell(metaboliteNameList))
    tmp = metaboliteNameList;
    clear metaboliteNameList;
    metaboliteNameList{1} = tmp;
end

for i = 1:length(metaboliteNameList)
    if nargin == 2
        rxnName = ['DM_' metaboliteNameList{i}];
        rxnNames{i}=rxnName;
    else
        rxnName = RxnName;
    end;
    metaboliteList = {metaboliteNameList{i}};
    model = addReaction(model,rxnName,metaboliteList,-1,false,0,1000,0,'Demand');
end
