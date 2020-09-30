function [solutionDel,D_square] = MOMA(modelWT,modelDel,solutionWTx,omega)

% Calculates MOMA solution to given wildtype solution

% modelWT:              Wildtype model
% modelDel:             Mutant model
% solutionWTx:          Wildtype solution (.x vector)
% omega (OPTIONAL):     weight vector (ones if omitted)


RxnsCheck = isequal(modelWT.rxns,modelDel.rxns);
MetsCheck = isequal(modelWT.mets,modelDel.mets);

if ~RxnsCheck || ~MetsCheck
    warning('WT and KO models must have identical reactions & metabolites');
end;



M = length(modelDel.mets);
N = length(modelDel.rxns);

if nargin == 3
    omega = ones(N,1);
end;

b = zeros(M,1);
A = modelDel.S;
c = - omega .* solutionWTx;
Q = diag(omega);
lb = modelDel.lb;
ub = modelDel.ub;
csense(1:M) = 'E';


[QPproblem.A,QPproblem.b,QPproblem.F,QPproblem.c,QPproblem.lb,QPproblem.ub,...
    QPproblem.csense,QPproblem.osense] = deal(A,b,Q,c,lb,ub,csense,1);

QPsolution = solveCobraQP(QPproblem);


if (QPsolution.stat > 0)
    solutionDel.x = QPsolution.full;
    D_square = sum((solutionWTx-solutionDel.x).^2);
end
solutionDel.stat = QPsolution.stat;
solutionDel.solver = QPsolution.solver;
solutionDel.time = QPsolution.time;


end