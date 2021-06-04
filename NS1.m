
%% This is the (NS-II) formulations in the paper 
%
% Optimal Network Slicing for Service-Oriented
% Networks with Flexible Routing and Guaranteed E2E Latency 
%
% by Wei-Kun Chen, Ya-Feng Liu, Antonio De Domenico, Zhi-Quan Luo, and Yu-Hong Dai
% to appear in IEEE Transactions on Network and Service Management, 2021 

% Send your comments and suggestions to   
%             chenweikun@bit.edu.cn        


function [sol,x,y,r,z,objvalue,sumy,Aver_nodedelay,Aver_linkdelay,Sollinkdelay,...
    Solnodedelay,Aver_path_num,Aver_max_path_num,totalr,Aver_path_ratio,Aver_min_path_ratio] = ...
    NS1(nP,nK1,isdelaycons,isdelayobj,originalgraphfile,virtualgraphfile,flowfile)

load(originalgraphfile);
load(virtualgraphfile);
load(flowfile);

%Parameter
nS = 3;
nK = nK1;
Profun=2;
newnV = nS*nV;
newnL = nL + 2*nS*nV;
newnI = nS*nV+nI;

%Placement varibales
y=binvar(nV,1);
x=binvar(newnV,nK,nS,'full');
x0=binvar(nV,nK,nS,'full');

%Routing variables
lambda=sdpvar(nK,nS+1,nP,'full');
r=sdpvar(newnL,nK,nS+1,nP,'full');
z=binvar(newnL,nK,nS+1,nP,'full');
theta=sdpvar(nK,nS+1,'full');
theta_L=sdpvar(nK,1);
theta_N=sdpvar(nK,1);

% Additional linearized variables
w1=sdpvar(newnV,nK,nS,nP,'full');
w2=sdpvar(newnV,nK,nS,nP,'full');

%bounds constraint
variablecons=[lambda>=0,r>=0,theta>=0,w1>=0,w2>=0,theta_L>=0,theta_N>=0];


%% VNF placement constraints
%Set x0[v,k,s] to be 0 if cloud node v cannot process function f_s^k
%We require that node V[1] can process all service functions
noprovidecons=[];
for k=1:nK
    for s=1:nS
        for v=2:nV
            curprocessfunctions=V_processfunction(1:Profun,v);
            if sum(FunctionChainC(k,s)==curprocessfunctions) == 0
                noprovidecons=noprovidecons+[x0(v,k,s)==0];
            end
        end
    end
end


%We require that at most one function in an SFC can be processed by the same node
atmostonecons=[];
for k=1:nK
    for v=1:newnV
        addition=0;
        for s=1:nS
            addition=addition+x(v,k,s);
        end
        atmostonecons=atmostonecons+[addition<=1];
    end
end

% Relate x and x0
xx0relationcons=[];
for k=1:nK
    for s=1:nS
        for v=1:nV
            addition=0;
            for t=1:nS
                addition=addition+x((v-1)*nS+t,k,s);
            end
            xx0relationcons=xx0relationcons+[x0(v,k,s)==addition];
        end
    end
end

%Every function should be processed by exactly one cloud node
everyfuncons=[];
for k=1:nK
    for s=1:nS
        addition=0;
        for v=1:nV
            addition=addition+x0(v,k,s);
        end
        everyfuncons=everyfuncons+[addition==1];
    end
end

% Relate x0 and y
x0yrelationcons=[];
for k=1:nK
    for s=1:nS
        for v=1:nV
            x0yrelationcons=x0yrelationcons+[x0(v,k,s)<=y(v)];
        end
    end
end

%Node capacity constraint
nodecapcons=[];
for v=1:nV
    addition=0;
    for k=1:nK
        for s=1:nS
            addition=addition+DataRate(k)*x0(v,k,s);
        end
    end
    nodecapcons=nodecapcons+[addition<=Node_cap(v)*y(v)];
end

%Total capacity constraint
totalcapcons=[];
addition=0;
for v=1:nV
    addition=addition+Node_cap(v)*y(v);
end
TotalDataRate=sum(DataRate(1:nK))*nS;
totalcapcons=totalcapcons+[addition>=TotalDataRate];


%% Routing constraints
%Enforce that the total data rates between two nodes hosting functions f_s^k and f_{s+1}^k to be equal to DataRate(k)
dataratecons=[];
for k=1:nK
    for s=1:(nS+1)
        addition=0;
        for p=1:nP
            addition=addition+lambda(k,s,p);
        end
        dataratecons=dataratecons+[addition==DataRate(k)];
    end
end

%at most at most one of the two links associated with (virtual) cloud node v is used by the p-th path of flow (k,s)
atmostonelinkcons=[];
for l=(nL+1):2:newnL
    atmostonelinkcons=atmostonelinkcons+[z(l,:,:,:)+z(l+1,:,:,:)<=1];
end

%relate r and z
rzrelationcons=[];
for k=1:nK
    rzrelationcons=rzrelationcons+[r(:,k,:,:)<=DataRate(k)*z(:,k,:,:)];
end

%Link capacity constraint
linkcapcons=[];
for l=1:nL
    linkcapcons=linkcapcons+[sum(sum(sum(r(l,:,:,:))))<=Link_cap(l)];
end

%Linearized constraint
% w1(v,k,s,p) = lambda(k,s,p)*x(v,k,s)
% w2(v,k,s,p) = lambda(k,s+1,p)*x(v,k,s)
linearizedcons=[];
for k=1:nK
    for s=1:nS
        for p=1:nP
            for v=1:newnV
                linearizedcons=linearizedcons+[w1(v,k,s,p)<=DataRate(k)*x(v,k,s)];
                linearizedcons=linearizedcons+[w2(v,k,s,p)<=DataRate(k)*x(v,k,s)];
                linearizedcons=linearizedcons+[w1(v,k,s,p)<=lambda(k,s,p)];
                linearizedcons=linearizedcons+[w2(v,k,s,p)<=lambda(k,s+1,p)];
                linearizedcons=linearizedcons+[w1(v,k,s,p)>=DataRate(k)*x(v,k,s)+lambda(k,s,p)-DataRate(k)];
                linearizedcons=linearizedcons+[w2(v,k,s,p)>=DataRate(k)*x(v,k,s)+lambda(k,s+1,p)-DataRate(k)];
            end
        end
    end
end

% SFC constraints
SFCcons2=[];
for k=1:nK
    for s=1:(nS+1)
        for i=1:newnI
            outgoinglinks=find(newL(:,1)==i);
            incominglinks=find(newL(:,2)==i);
            addition1=zeros(1,1,1,nP);
            addition2=zeros(1,1,1,nP);
            for t=1:length(outgoinglinks)
                addition1=addition1+r(outgoinglinks(t),k,s,:);
                addition2=addition2+z(outgoinglinks(t),k,s,:);
            end
            addition3=zeros(1,1,1,nP);
            addition4=zeros(1,1,1,nP);
            for t=1:length(incominglinks)
                addition3=addition3+r(incominglinks(t),k,s,:);
                addition4=addition4+z(incominglinks(t),k,s,:);
            end
            
            if s==1
                if i==PairC(k,1)
                    SFCcons2=SFCcons2+[addition1==lambda(k,s,:)];
                    SFCcons2=SFCcons2+[addition2==ones(1,1,1,nP)];
                    SFCcons2=SFCcons2+[addition3==zeros(1,1,1,nP)];
                    SFCcons2=SFCcons2+[addition4==zeros(1,1,1,nP)];
                elseif i>=nI+1
                    v = i-nI;
                    SFCcons2=SFCcons2+[addition1==zeros(1,1,1,nP)];
                    SFCcons2=SFCcons2+[addition2==zeros(1,1,1,nP)];
                    SFCcons2=SFCcons2+[addition3==ones(1,1,1,nP).*w1(v,k,s,:)];
                    SFCcons2=SFCcons2+[addition4==ones(1,1,1,nP)*x(v,k,s)];
                else
                    SFCcons2=SFCcons2+[addition1==addition3];
                    SFCcons2=SFCcons2+[addition2==addition4];
                end
            elseif s==(nS+1)
                if i==PairC(k,2)
                    SFCcons2=SFCcons2+[addition1==zeros(1,1,1,nP)];
                    SFCcons2=SFCcons2+[addition2==zeros(1,1,1,nP)];
                    SFCcons2=SFCcons2+[addition3==lambda(k,s,:)];
                    SFCcons2=SFCcons2+[addition4==ones(1,1,1,nP)];
                elseif i>=nI+1
                    v=i-nI;
                    SFCcons2=SFCcons2+[addition1==ones(1,1,1,nP).*w2(v,k,s-1,:)];
                    SFCcons2=SFCcons2+[addition2==ones(1,1,1,nP)*x(v,k,s-1)];
                    SFCcons2=SFCcons2+[addition3==zeros(1,1,1,nP)];
                    SFCcons2=SFCcons2+[addition4==zeros(1,1,1,nP)];
                else
                    SFCcons2=SFCcons2+[addition1==addition3];
                    SFCcons2=SFCcons2+[addition2==addition4];
                end
            else
                if i>=nI+1
                    v=i-nI;
                    SFCcons2=SFCcons2+[addition1==ones(1,1,1,nP).*w2(v,k,s-1,:)];
                    SFCcons2=SFCcons2+[addition2==ones(1,1,1,nP)*x(v,k,s-1)];
                    SFCcons2=SFCcons2+[addition3==ones(1,1,1,nP).*w1(v,k,s,:)];
                    SFCcons2=SFCcons2+[addition4==ones(1,1,1,nP)*x(v,k,s)];
                else
                    SFCcons2=SFCcons2+[addition1==addition3];
                    SFCcons2=SFCcons2+[addition2==addition4];
                end
            end
        end
    end
end

%% Delay constraint
%Function-function delay
mediadelaycons=[];
for k=1:nK
    for s=1:(nS+1)
        for p=1:nP
            addition=0;
            for l=1:nL
                addition=addition+Link_delay(l)*z(l,k,s,p);
            end
            mediadelaycons=mediadelaycons+[addition<=theta(k,s)];
        end
    end
end

%Link delay constraint
linkdelaycons=[];
for k=1:nK
    linkdelaycons=linkdelaycons+[sum(theta(k,:))==theta_L(k)];
end

%Node delay constraint
nodedelaycons=[];
for k=1:nK
    addition=0;
    for v=1:nV
        addition=addition+sum(x0(v,k,:)*Node_delay(v));
    end
    nodedelaycons=nodedelaycons+[addition==theta_N(k)];
end

%E2E delay constraint
E2Econs=[];
if isdelaycons
    for k=1:nK
        E2Econs=E2Econs+[theta_L(k)+theta_N(k)<=E2Edelays(k)];
    end
end

% Objective function
if isdelayobj == 1
    if isdelaycons == 1
        obj = sum(y) + 0.001*(sum(theta_L)+sum(theta_N));
    else
        addition=0;
        for i=1:nL
            addition=addition+sum(sum(sum(r(l,:,:,:))));
        end
        %       nL
        obj = sum(y) + 0.0001*sum(sum(sum(sum(r(1:nL,:,:,:)))));
    end
else
    obj = sum(y);
end


ops = sdpsettings('solver','gurobi');
ops.gurobi.TimeLimit = 600;
ops.gurobi.MIPGap = 0.000;

sol = optimize([variablecons, noprovidecons, atmostonecons, xx0relationcons,...
    everyfuncons,x0yrelationcons,nodecapcons,totalcapcons,dataratecons,...
    atmostonelinkcons, rzrelationcons, linkcapcons,...
    linearizedcons,SFCcons2,...
    mediadelaycons,linkdelaycons,nodedelaycons,E2Econs],obj,ops);
objvalue=0;
valuey=value(y);
sumy=sum(valuey);
Aver_nodedelay=0;
Aver_linkdelay=0;
Aver_path_num=0;
Aver_max_path_num=0;
Sollinkdelay=0;
Solnodedelay=0;
totalr=0;
Aver_path_ratio=0;
Aver_min_path_ratio = 0;


if sumy~=0
    objvalue=value(obj);
    Solnodedelay=value(theta_N);
    valuer=value(r);
    valuez=value(z);
    valuelambda=value(lambda);
    totalr = sum(sum(sum(sum(valuer(1:nL,:,:,:)))));
    
    %value(sum(theta_N))
    %value(x0)
    %value(theta_L)
    
    
    %compute real delay
    Sollinkdelay=zeros(nK,1);
    for k=1:nK
        for s=1:(nS+1)
            tmp=0;
            for p=1:nP
                tmp=max(tmp,sum((valuer(:,k,s,p)>1e-6).*newLink_delay));
            end
            Sollinkdelay(k) = Sollinkdelay(k)+tmp;
        end
    end
    
    %check whether this solution satisfies the delay constraints
    if isdelaycons==0
        Delay=zeros(nK,1);
        for k=1:nK
            Delay(k)=Sollinkdelay(k)+Solnodedelay(k);
            if (Delay(k)>E2Edelays(k))
                sumy=0;
                break;
            end
        end
    end
    % Average delays
    Aver_nodedelay=mean(Solnodedelay);
    Aver_linkdelay=mean(Sollinkdelay);
    
    %compute the number of paths and the ratio of data rate on each path
    Aver_path_num=1;
    Aver_max_path_num=1;
    if isdelaycons==1 && nP==2
        eps=1.0e-6;
        path=zeros(nK,1);
        minrateonpath=zeros(nK,1);
        for k=1:nK
            R=zeros(nL,1,1,1);
            AllR=zeros(2*(nS+1),1);
            consL=zeros(2*(nS+1),2);
            for s=1:(nS+1)
                consL(2*(s-1)+1,1)=s;
                consL(2*s,1)=s;
                consL(2*(s-1)+1,2)=s+1;
                consL(2*s,2)=s+1;
                % test the data rate on each path
                if valuelambda(k,s,1)<eps || valuelambda(k,s,2)<eps
                    AllR(2*(s-1)+1)=DataRate(k);
                    AllR(2*s)=0;
                else
                    tmp = abs(valuez(1:newnL,k,s,1)-valuez(1:newnL,k,s,2));
                    if sum(tmp) < eps
                        % the data rate on both paths are positive but indeed they are the same path
                        AllR(2*(s-1)+1)=DataRate(k);
                        AllR(2*s)=0;
                    else
                        AllR(2*(s-1)+1)=valuelambda(k,s,1);
                        AllR(2*s)=valuelambda(k,s,2);
                    end
                end
            end
            [path(k) minrateonpath(k)]=paths_rates(DataRate(k),nS+2,2*(nS+1),consL,AllR,1,5);
        end
        Aver_path_num=mean(path);
        Aver_max_path_num=max(path);
        Aver_path_ratio = mean(minrateonpath);
        Aver_min_path_ratio = min(minrateonpath);
    end
else
    
end

end
