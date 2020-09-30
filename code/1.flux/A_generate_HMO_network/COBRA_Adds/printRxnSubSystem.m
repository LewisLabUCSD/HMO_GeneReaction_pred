function subSysRxnList = printRxnSubSystem(model,subSystem,printFlag)

if nargin == 2
    printFlag = true;
end;

if ~isempty(subSystem)
    subSysSel = find(strcmp(model.subSystems,subSystem));
else
    subSysSel = cellfun(@isempty,model.subSystems);
end

subSysRxnList = model.rxns(subSysSel);
if printFlag
    printRxnFormula(model,subSysRxnList);
end;


end