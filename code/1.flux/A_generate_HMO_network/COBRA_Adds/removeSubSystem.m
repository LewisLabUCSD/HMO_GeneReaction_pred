function newmodel = removeSubSystem(model,subSystem)

if ~isempty(subSystem)
    subSysSel = strcmp(model.subSystems,subSystem);
else
    subSysSel = cellfun(@isempty,model.subSystems);
end

subSysRxnList = model.rxns(subSysSel);
newmodel = removeRxns(model,subSysRxnList);


end