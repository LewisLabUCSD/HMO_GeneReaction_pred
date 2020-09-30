% Code to generate all the pre-combination of reactions that allow speed up
% the enumeration of sub-networks

load('Network_Red_ext');
    model=Network_Red;

    DFLNT = [190 191 192];
    FLNH = [194 195 196 197 198 199 200];
    DFLNH =[202 203 204 205 206 207 208];
    FDSLNH =[210 211 212 213 214 215 216];
    DSLNH=[218 219 220];
   
    combi=allcomb(DFLNT,FLNH,DFLNH,FDSLNH,DSLNH);
    
    All_rxns=[DFLNT,FLNH,DFLNH,FDSLNH,DSLNH];
    global_combi={};
    for i=1:size(combi,1)

        all_sol=[];

        i
        [neighborRxns] = findUpstreamRxns(model,model.rxns(combi(i,1)),'true');
        DFLNT_up=findRxnIDs(model,neighborRxns);
        [neighborRxns] = findUpstreamRxns(model,model.rxns(combi(i,2)),'true');
        FLNH_up=findRxnIDs(model,neighborRxns);
        [neighborRxns] = findUpstreamRxns(model,model.rxns(combi(i,3)),'true');
        DFLNH_up=findRxnIDs(model,neighborRxns);
        [neighborRxns] = findUpstreamRxns(model,model.rxns(combi(i,4)),'true');
        FDSLNH_up=findRxnIDs(model,neighborRxns);
        [neighborRxns] = findUpstreamRxns(model,model.rxns(combi(i,5)),'true');
        DSLNH_up=findRxnIDs(model,neighborRxns);
        
        combi_upstream = allcomb(DFLNT_up,FLNH_up,DFLNH_up,FDSLNH_up,DSLNH_up);

        
        for j=1:size(combi_upstream,1)
            [neighborRxns] = findUpstreamRxns(model,model.rxns(combi_upstream(j,1)),'true');
            Up_1=findRxnIDs(model,neighborRxns);
            [neighborRxns] = findUpstreamRxns(model,model.rxns(combi_upstream(j,2)),'true');
            Up_2=findRxnIDs(model,neighborRxns);
            [neighborRxns] = findUpstreamRxns(model,model.rxns(combi_upstream(j,3)),'true');
            Up_3=findRxnIDs(model,neighborRxns);
            [neighborRxns] = findUpstreamRxns(model,model.rxns(combi_upstream(j,4)),'true');
            Up_4=findRxnIDs(model,neighborRxns);
            [neighborRxns] = findUpstreamRxns(model,model.rxns(combi_upstream(j,5)),'true');
            Up_5=findRxnIDs(model,neighborRxns);

            Up_upstream = allcomb(Up_1,Up_2,Up_3,Up_4,Up_5);

            for k=1:size(Up_upstream,1)
                all_sol(end+1,:)=[combi(i,:) combi_upstream(j,:) Up_upstream(k,:)];
            end
        end
    global_combi{i}=all_sol;
    end
   