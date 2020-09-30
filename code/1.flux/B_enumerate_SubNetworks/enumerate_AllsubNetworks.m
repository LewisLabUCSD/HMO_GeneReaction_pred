function [ SOL, perf] = enumerate_AllsubNetworks(p)
    %UNTITLED Summary of this function goes here
    %   Detailed explanation goes here
   % changeCobraSolver('gurobi5', 'LP');
   % changeCobraSolver('gurobi5', 'MILP');
    
    changeCobraSolver('glpk', 'LP');
    changeCobraSolver('glpk', 'MILP');
    load('Network_Red_ext');
    model=Network_Red;
    load('global_combi');
    
    DFLNT = [190 191 192];
    FLNH = [194 195 196 197 198 199 200];
    DFLNH =[202 203 204 205 206 207 208];
    FDSLNH =[210 211 212 213 214 215 216];
    DSLNH=[218 219 220];
   
    %combi=allcomb(DFLNT,FLNH,DFLNH,FDSLNH,DSLNH);
    
    All_rxns=[DFLNT,FLNH,DFLNH,FDSLNH,DSLNH];

    Total=0;
    objective=[];
    [m,n] = size(model.S);
    impossible=[];
    SOL=[];
    

        model_init=Network_Red;
        model_init.lb(All_rxns)=0;
        model_init.ub(All_rxns)=0;

        combi=global_combi{p};


        for C_up=1:size(combi,1)
            model=model_init;
            model.lb(combi(C_up,:))=1e-3;
            model.ub(combi(C_up,1:5))=10;


                
            maxFlux=zeros(1,221);
            minFlux=zeros(1,221);
            
            for k=1:221
            	model.c=zeros(221,1);
            	model.c(k)=1;
            	FBAsol = optimizeCbModel(model,'max');
                maxFlux(k)=FBAsol.f;
            end
            model.ub(find(maxFlux==0))=0;

            % variables:
            %    v's (n), y's (n)

            % constraints:
            %    m mass balance constraints
            A = [model.S zeros(m,n)];
            % constrain UB fluxes w/ integer constraints
            A = [A; 
                [eye(n,n) -diag(model.ub)] ];
            % constrain LB fluxes w/ integer constraints
            A = [A; 
                eye(n,n) -diag(model.lb) ];

            % create structure
            MILPproblem.A = A;
            MILPproblem.b  = [zeros(m,1);zeros(2*n,1)];
            MILPproblem.c  = [zeros(n,1); ones(n,1)];
            MILPproblem.lb = [model.lb; zeros(n,1)];
            MILPproblem.ub = [model.ub; ones(n,1)];    
            MILPproblem.csense(1:m) ='E';
            MILPproblem.csense(m+1:m+n) ='L';
            MILPproblem.csense(m+n+1:m+2*n) ='G';

            %osense Objective sense (-1 max, +1 min);
            MILPproblem.osense = 1;
            MILPproblem.vartype(1:n) = 'C';
            MILPproblem.vartype(n+1:2*n) = 'I';
            MILPproblem.x0 = [];%zeros(2*n,1);
            MILPsol = solveCobraMILP(MILPproblem);

            maxObjective = MILPsol.obj;

            objective=[objective; maxObjective];

            if maxObjective==0
                %C_up=C_up+1;
                continue
            end

            prevNZ = abs(MILPsol.cont) > .0001;
            NZ = prevNZ;

            solution.fluxes=[];
            solution.nonzero =[];

            MILPsol2=MILPsol;
            NZ = abs(MILPsol.cont) > .0001;

            prevNZ=  NZ(:,end) ;
            NZ(:,end+1) = abs(MILPsol2.full(1:n))>.000000001;
            PrevNZ = NZ(:,end);

            count=1;
            while count==1
                MILPproblem2=MILPproblem;

                   % variables:
                %    v's (n), y's (n) w's (n)  3n total variables

                    % constraints:
                %   resize of m mass balance constraints for w's addition
                MILPproblem2.A = [MILPproblem2.A zeros(m+2*n,n)];


                % constrain w+y <=1 
                MILPproblem2.A = [MILPproblem2.A; zeros(n,n), eye(n,n), eye(n,n) ];
                MILPproblem2.b = [MILPproblem2.b; ones(n,1)];
                MILPproblem2.csense(m+2*n+1:m+3*n) = 'L';


                % constrain with previous zero results
                %MILPproblem2.A = [MILPproblem2.A; zeros(1,n), prevNZ', zeros(1,n) ];
                MILPproblem2.A = [MILPproblem2.A; zeros(1,n) zeros(1,n) prevNZ' ];
                MILPproblem2.b = [MILPproblem2.b; 1];
                MILPproblem2.csense(m+(3*n)+1:m+(3*n)+1) = 'G';

                % constrain with previous results (altbases)
                for nI = 1:size(NZ,2)
                    %MILPproblem2.A = [MILPproblem2.A; [zeros(1,n), zeros(1,n) NZ(:,i)']];
                    MILPproblem2.A = [MILPproblem2.A; [zeros(1,n) NZ(:,nI)', zeros(1,n) ]];
                    MILPproblem2.b(end+1) = sum(NZ(:,nI))-1;
                    MILPproblem2.csense(end+1) = 'L';
                end


                % vartype
                for nJ = 1:n
                    MILPproblem2.vartype(end+1) = 'B';
                end

                % lb,ub
                MILPproblem2.lb = [MILPproblem.lb; zeros(n,1)];
                MILPproblem2.ub = [MILPproblem.ub; ones(n,1)];    
                % c
                MILPproblem2.c = [MILPproblem.c; zeros(n,1)];

                MILPsol2 = solveCobraMILP(MILPproblem2);
                if (abs(MILPsol2.obj)==0)
                	count=2;
                    continue
                end
                
                NZ(:,end+1) = abs(MILPsol2.full(1:n))>.000000001;
                PrevNZ = NZ(:,end);

                solution.fluxes = [solution.fluxes,MILPsol2.full(1:n)];
                solution.nonzero = [solution.nonzero, NZ(:,end)];
  
                Total=Total+1;
            end
                SOL =[SOL solution.nonzero];
                display(['Combi tested = ',num2str(p),' and upstream combinations tested = ',num2str(C_up),'/',num2str(size(combi,1))])
                perf=[num2str(C_up),'/',num2str(size(combi,1))];
        end

end


