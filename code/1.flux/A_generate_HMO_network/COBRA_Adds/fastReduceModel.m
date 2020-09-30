function [modelRed,mins,maxes,hasFlux] = fastReduceModel(model,tol,WCount,changeBoundsFlag,checkConsistencyFlag)
% Model reduction using fastFVA

    function modelOK = checkConsistency(model,modelRed,tol)
    % modelOK = checkConsistency(model,modelRed,tol)
    if (sum(model.c ~= 0) > 0)
        % Original model
        solOrigMax = optimizeCbModel(model,'max',[]);
        solOrigMin = optimizeCbModel(model,'min',[]);
        % Reduced model
        solRedMax = optimizeCbModel(modelRed,'max',[]);
        solRedMin = optimizeCbModel(modelRed,'min',[]);
        diffMax = abs(solRedMax.f - solOrigMax.f);
        diffMin = abs(solRedMin.f - solOrigMin.f);
        if (diffMax > tol || diffMin > tol)
            fprintf('Inconsistent objective values %g %g %g %g\n',solOrigMax.f,solRedMax.f,solOrigMin.f,solRedMin.f);
            modelOK = false;
        else
            fprintf('Model is consistent\n');
            modelOK = true;
        end
    else
        modelOK = true;
    end
    end


% ------------------------------------------------------------------------

if WCount > 1
    SetWorkerCount(WCount);
end
[mins,maxes] = fastFVA(model,100,'max','glpk');
%[mins,maxes] = fastFVA(model,100,'max','cplex');

%fastFVA does not work

mins = mins .* (abs(mins)>tol);
maxes = maxes .* (abs(maxes)>tol);

% Create a list of flux indexes that have non-zero flux (hasFlux)
hasFluxSel = (maxes ~= 0 | mins ~= 0);
hasFlux = find(hasFluxSel);
hasFlux = columnVector(hasFlux);


% Reduce model
modelRed = removeRxns(model,model.rxns(~hasFluxSel));

% Update bounds
if (changeBoundsFlag)
    modelRed.lb = columnVector(mins(hasFlux));
    modelRed.ub = columnVector(maxes(hasFlux));
    selInconsistentBounds = (modelRed.ub < modelRed.lb);
    modelRed.ub(selInconsistentBounds) = modelRed.lb(selInconsistentBounds);
end

% Consistency check
checkConsistencyFlag = false; % DOES NOT MAKE SENSE IF changeBoundsFlag = true
if (checkConsistencyFlag)
    fprintf('Perform model consistency check\n');
    checkConsistency(model,modelRed,tol);
end

% if WCount > 1
%     SetWorkerCount(1)
% end;

end