function MetSinkIDs = findSinkIDs(MetList,model,subSys)
% Find sink reaction IDs (located in subSys) in model for metabolites listed in MetList

if ~iscell(MetList)
    MetList = {MetList};
end;
MetIDs = findMetIDs(model,MetList);
subSysIDs = find(strcmp(model.subSystems,subSys));
L = length(MetList);
MetSinkIDs = zeros(1,length(MetIDs));

for i = 1:L
    MetSel = find(model.S(MetIDs(i),subSysIDs));
    MetSinkIDs(i) = subSysIDs(MetSel);
end;

end