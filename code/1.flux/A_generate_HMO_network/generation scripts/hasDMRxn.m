function [YesOrNo,DMRxnID] = hasDMRxn(model,met)

RxnList = findRxnsFromMets(model,met);
RxnIDs = findRxnIDs(model,RxnList);
YesOrNo = false;
for j = 1:length(RxnIDs)
    if length(find(model.S(:,RxnIDs(j)))) == 1
        YesOrNo = true;
        DMRxnID = RxnIDs(j);
    end;
end;

end
