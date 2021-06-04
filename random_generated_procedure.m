% This file descirbes the procedure to generate the random networks in the
% paper
%
% Optimal Network Slicing for Service-Oriented
% Networks with Flexible Routing and Guaranteed E2E Latency 
%
% by Wei-Kun Chen, Ya-Feng Liu, Antonio De Domenico, Zhi-Quan Luo, and Yu-Hong Dai
% to appear in IEEE Transactions on Network and Service Management, 2021 

% Send your comments and suggestions to   
%             chenweikun@bit.edu.cn        


load('r_seed.mat');
rng(r_seed);

% number of vertices
nI=6;
I = 1:nI;

% number of cloud nodes
nV=3;

% length in the SFC
nS=3;

%sparsity ratio
ratio = 0.6;

%directory
random_network=['random_networks'];
%number of problems
nprob = 100;
%number of function types
nflow_types =5;
% number of flows
nK = 5;
%generate the combinations of the service function chains
flow_combs=combntns([1:nflow_types],nS);
nflow_combs=combntns(nflow_types,nS);


V_processfunction=zeros(nflow_types,nV);
DG=sparse(zeros(nI, nI));

for t=1:nprob
    % generate a strongly connected graph
    while 1
        x = rand(nI,1)*100;
        y = rand(nI,1)*100;
        Link = rand(nI,nI) > (1-ratio);
        D = round(squareform(pdist([x y])));
        
        for i=1:nI
            Link(i,i) = 0;
        end
        
        nL=sum(sum(Link));
        Edge_cap = 0.5 + 3*rand(nI,nI);
        DG = sparse(D.*Link);
        % strongly connected
        [S,C] = graphconncomp(DG);
        
        if S==1
            break;
        end
    end
    
    %links
    L=zeros(nL,2);
    %edge_cap
    Link_cap=zeros(nL,1);
    %edge_delay
    Link_delay=zeros(nL,1);
    
    dist=graphallshortestpaths(DG);
    aver_shortestpaths=sum(sum(dist))/(nI*(nI-1));
    
    DG=DG/aver_shortestpaths;
    is_V=(rand(1,nI)<=0.5);
    while sum(is_V)~=nV
        is_V=(rand(1,nI)<0.5);
    end
    
    V = I(is_V);
    
    originalgraphfile=[random_network '/graph' '-' num2str(t)];
    
    k=1;
    for i=1:nI
        for j=1:nI
            if(Link(i,j)==1)
                L(k,1) = i;
                L(k,2) = j;
                Link_delay(k)=D(i,j)/aver_shortestpaths;
                Link_cap(k)= Edge_cap(i,j);
                k=k+1;
            end
        end
    end
    
    Node_cap=round(6+6*rand(1,nV));
    Node_delay=0.8+0.4*rand(1,nV);
    for v=1:nV
        tmp1=1:nflow_types;
        [tmp,index]=sort(rand(1,nflow_types));
        V_processfunction(:,v) = tmp1(index);
    end
    % save the graph data
    save(originalgraphfile,'nI','nL','nV','V','V_processfunction','L','Link_delay','Link_cap','Node_delay','Node_cap');
    
    FunctionChainC = zeros(nK,nS);
    PairC = zeros(nK,2);
    
    % generate source-destinction paris
    nnonV=nI-nV;
    nonV=I-I.*is_V;
    res_nonV=nonV(nonV>0);
    nonNFV_combs=combntns(res_nonV,2);
    tmp = [nonNFV_combs(:,2) nonNFV_combs(:,1)];
    nonNFV_combs=[nonNFV_combs;tmp];
    nnonNFV_combs=2*combntns(nnonV,2);
    tmp=round(0.5+(nnonNFV_combs-0.5)*rand(nK,1));
    PairC=nonNFV_combs(tmp',:);
    
    % generate service function chains
    tmp=round(0.5+(nflow_combs-0.5)*rand(nK,1));
    flowsfunctions=flow_combs(tmp',:);
    for j=1:nK
        d=randperm(nS);
        tmp=flowsfunctions(j,:);
        FunctionChainC(j,:)=tmp(d);
    end
    
    
    DataRate=ones(nK,1);
    shortpathdist=zeros(nK,1);
    for i=1:nK
        [shortpathdist(i),path,pred]=graphshortestpath(DG,PairC(i,1),PairC(i,2),'Directed','true');
    end
    E2Edelays = round(3+(6*shortpathdist)+2*rand());
    
    flowfile=[random_network '/SFCs' '-' num2str(t)];
    
    save(flowfile,'nK','nS','nflow_types','PairC','E2Edelays','DataRate','FunctionChainC');
    
    filename=[random_network '/virtualgraph' '-' num2str(t)];
    constructed_virtualgraph(nI,nS,nV,nL,V,L,Link_cap,Link_delay,filename);
    
end
















