
%% This is the (NS-I) formulations in the paper 
%
% Optimal Network Slicing for Service-Oriented
% Networks with Flexible Routing and Guaranteed E2E Latency 
%
% by Wei-Kun Chen, Ya-Feng Liu, Antonio De Domenico, Zhi-Quan Luo, and Yu-Hong Dai
% to appear in IEEE Transactions on Network and Service Management, 2021 

% Send your comments and suggestions to   
%             chenweikun@bit.edu.cn        


function [sol,objvalue] = NS_I(nP,nK1,isdelaycons,isdelayobj,...
    originalgraphfile,virtualgraphfile,flowfile)


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

%hosting pair
npaircloudnodes=(newnV)*(newnV-1);
paircloudnodes=zeros(npaircloudnodes,2);

w=binvar(npaircloudnodes,nK,nS-1,'full');

q=1;
for i=1:newnV
    for j=1:newnV
        if i==j
            continue;
        else
            paircloudnodes(q,1)=i;
            paircloudnodes(q,2)=j;
            q=q+1;
        end
    end
end

%routing variables
R=sdpvar(nK,nS-1,npaircloudnodes,nP,'full');
z=binvar(newnL,nK,nS-1,npaircloudnodes,nP,'full');
r=sdpvar(newnL,nK,nS-1,npaircloudnodes,nP,'full');

R_o=sdpvar(newnV,nK,nP,'full');
z_o=binvar(newnL,newnV,nK,nP,'full');
r_o=sdpvar(newnL,newnV,nK,nP,'full');

R_d=sdpvar(newnV,nK,nP,'full');
z_d=binvar(newnL,newnV,nK,nP,'full');
r_d=sdpvar(newnL,newnV,nK,nP,'full');

theta=sdpvar(nK,nS+1,'full');
theta_L=sdpvar(nK,1);
theta_N=sdpvar(nK,1);


%bounds constraint
variablecons=[R>=0,R_o>=0,R_d>=0,r>=0,r_o>=0,r_d>=0,...
    theta>=0,theta_L>=0,theta_N>=0];


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


%Linearized constraint
% w(v_s,v_{s+1},k,s) = x(v_s,k,s)*x(v_{s+1},k,s)
linearizedcons=[];
for k=1:nK
    for s=1:(nS-1)
        for q=1:npaircloudnodes
            linearizedcons=linearizedcons+[w(q,k,s)<=x(paircloudnodes(q,1),k,s)];
            linearizedcons=linearizedcons+[w(q,k,s)<=x(paircloudnodes(q,2),k,s+1)];
            linearizedcons=linearizedcons+[w(q,k,s)>=x(paircloudnodes(q,1),k,s)+...
                x(paircloudnodes(q,2),k,s+1)-1];
        end
    end
end


%%Enforce that the total data rates between two nodes hosting functions f_s^k and f_{s+1}^k to be equal to DataRate(k)
dataratecons=[];

%O(k)--> f^1 %f^{nS}->D(k)
for k=1:nK
    for v=1:newnV
        addition1=0;
        addition2=0;
        for p=1:nP
            addition1=addition1+R_o(v,k,p);
            addition2=addition2+R_d(v,k,p);
        end
        dataratecons=dataratecons+[addition1==DataRate(k)*x(v,k,1)];
        dataratecons=dataratecons+[addition2==DataRate(k)*x(v,k,nS)];
    end
end

%f^s->f^{s+1}
for k=1:nK
    for s=1:(nS-1)
        for q=1:npaircloudnodes
            addition=0;
            for p=1:nP
                addition=addition+R(k,s,q,p);
            end
            dataratecons=dataratecons+[addition==DataRate(k)*w(q,k,s)];
        end
    end
end


%related x, r, and z
xrzrelationcons=[];
tmp=ones(newnL,1,1,nP);
for k=1:nK
    for v=1:newnV
        xrzrelationcons=xrzrelationcons+[z_o(:,v,k,:)<=x(v,k,1)*tmp];
        xrzrelationcons=xrzrelationcons+[r_o(:,v,k,:)<=z_o(:,v,k,:)*DataRate(k)];
        xrzrelationcons=xrzrelationcons+[z_d(:,v,k,:)<=x(v,k,nS)*tmp];
        xrzrelationcons=xrzrelationcons+[r_d(:,v,k,:)<=z_d(:,v,k,:)*DataRate(k)];
    end
end
tmp=ones(newnL,1,1,1,nP);
for k=1:nK
    for s=1:(nS-1)
        for q=1:npaircloudnodes
            xrzrelationcons=xrzrelationcons+[z(:,k,s,q,:)<=w(q,k,s)*tmp];
            xrzrelationcons=xrzrelationcons+[r(:,k,s,q,:)<=z(:,k,s,q,:)*DataRate(k)];
        end
    end
end


%Link capacity constraint
linkcapcons=[];
for l=1:nL
    addition=sum(sum(sum(r_o(l,:,:,:))))+sum(sum(sum(r_d(l,:,:,:))))+...
        sum(sum(sum(sum(r(l,:,:,:,:)))));
    linkcapcons=linkcapcons+[addition<=Link_cap(l)];
end

SFCcons2=[];
%source node
for k=1:nK
    for v=1:newnV
        for i=1:newnI
            outgoinglinks=find(newL(:,1)==i);
            incominglinks=find(newL(:,2)==i);
            addition1=zeros(1,1,1,nP);
            addition2=zeros(1,1,1,nP);
            for t=1:length(outgoinglinks)
                addition1=addition1+r_o(outgoinglinks(t),v,k,:);
                addition2=addition2+z_o(outgoinglinks(t),v,k,:);
            end
            addition3=zeros(1,1,1,nP);
            addition4=zeros(1,1,1,nP);
            for t=1:length(incominglinks)
                addition3=addition3+r_o(incominglinks(t),v,k,:);
                addition4=addition4+z_o(incominglinks(t),v,k,:);
            end
            if i==PairC(k,1)
                SFCcons2=SFCcons2+[addition1==R_o(v,k,:)];
                SFCcons2=SFCcons2+[addition2==ones(1,1,1,nP)*x(v,k,1)];
                SFCcons2=SFCcons2+[addition3==zeros(1,1,1,nP)];
                SFCcons2=SFCcons2+[addition4==zeros(1,1,1,nP)];
            elseif i==nI+v
                SFCcons2=SFCcons2+[addition1==zeros(1,1,1,nP)];
                SFCcons2=SFCcons2+[addition2==zeros(1,1,1,nP)];
                SFCcons2=SFCcons2+[addition3==R_o(v,k,:)];
                SFCcons2=SFCcons2+[addition4==ones(1,1,1,nP)*x(v,k,1)];
            else
                SFCcons2=SFCcons2+[addition1==addition3];
                SFCcons2=SFCcons2+[addition2==addition4];
            end
        end
    end
end

%destincnation node
for k=1:nK
    for v=1:newnV
        for i=1:newnI
            outgoinglinks=find(newL(:,1)==i);
            incominglinks=find(newL(:,2)==i);
            addition1=zeros(1,1,1,nP);
            addition2=zeros(1,1,1,nP);
            for t=1:length(outgoinglinks)
                addition1=addition1+r_d(outgoinglinks(t),v,k,:);
                addition2=addition2+z_d(outgoinglinks(t),v,k,:);
            end
            addition3=zeros(1,1,1,nP);
            addition4=zeros(1,1,1,nP);
            for t=1:length(incominglinks)
                addition3=addition3+r_d(incominglinks(t),v,k,:);
                addition4=addition4+z_d(incominglinks(t),v,k,:);
            end
            
            if i==nI+v
                SFCcons2=SFCcons2+[addition1==R_d(v,k,:)];
                SFCcons2=SFCcons2+[addition2==ones(1,1,1,nP)*x(v,k,nS)];
                SFCcons2=SFCcons2+[addition3==zeros(1,1,1,nP)];
                SFCcons2=SFCcons2+[addition4==zeros(1,1,1,nP)];
            elseif i==PairC(k,2)
                SFCcons2=SFCcons2+[addition1==zeros(1,1,1,nP)];
                SFCcons2=SFCcons2+[addition2==zeros(1,1,1,nP)];
                SFCcons2=SFCcons2+[addition3==R_d(v,k,:)];
                SFCcons2=SFCcons2+[addition4==ones(1,1,1,nP)*x(v,k,nS)];
            else
                SFCcons2=SFCcons2+[addition1==addition3];
                SFCcons2=SFCcons2+[addition2==addition4];
            end
        end
    end
end

%intermediate nodes
for k=1:nK
    for s=1:(nS-1)
        for q=1:npaircloudnodes
            for i=1:newnI
                outgoinglinks=find(newL(:,1)==i);
                incominglinks=find(newL(:,2)==i);
                addition1=zeros(1,1,1,1,nP);
                addition2=zeros(1,1,1,1,nP);
                for t=1:length(outgoinglinks)
                    addition1=addition1+r(outgoinglinks(t),k,s,q,:);
                    addition2=addition2+z(outgoinglinks(t),k,s,q,:);
                end
                addition3=zeros(1,1,1,nP);
                addition4=zeros(1,1,1,nP);
                for t=1:length(incominglinks)
                    addition3=addition3+r(incominglinks(t),k,s,q,:);
                    addition4=addition4+z(incominglinks(t),k,s,q,:);
                end
                v_s=paircloudnodes(q,1);
                v_s_1=paircloudnodes(q,2);
                
                if i==(nI+v_s)
                    SFCcons2=SFCcons2+[addition1==R(k,s,q,:)];
                    SFCcons2=SFCcons2+[addition2==ones(1,1,1,1,nP)*w(q,k,s)];
                    SFCcons2=SFCcons2+[addition3==zeros(1,1,1,1,nP)];
                    SFCcons2=SFCcons2+[addition4==zeros(1,1,1,1,nP)];
                elseif i==(nI+v_s_1)
                    SFCcons2=SFCcons2+[addition1==zeros(1,1,1,1,nP)];
                    SFCcons2=SFCcons2+[addition2==zeros(1,1,1,1,nP)];
                    SFCcons2=SFCcons2+[addition3==R(k,s,q,:)];
                    SFCcons2=SFCcons2+[addition4==ones(1,1,1,nP)*w(q,k,s)];
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
        if s==1
            for p=1:nP
                addition=0;
                for l=1:nL
                    addition=addition+Link_delay(l)*sum(z_o(l,:,k,p));
                end
                mediadelaycons=mediadelaycons+[addition<=theta(k,s)];
            end
        elseif s==(nS+1)
            for p=1:nP
                addition=0;
                for l=1:nL
                    addition=addition+Link_delay(l)*sum(z_d(l,:,k,p));
                end
                mediadelaycons=mediadelaycons+[addition<=theta(k,s)];
            end
        else
            for p=1:nP
                addition=0;
                for l=1:nL
                    addition=addition+Link_delay(l)*sum(z(l,k,s-1,:,p));
                end
                mediadelaycons=mediadelaycons+[addition<=theta(k,s)];
            end
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

%objective function
if isdelayobj == 1
    if isdelaycons == 1
        obj = sum(y) + 0.001*(sum(theta_L)+sum(theta_N));
    else
        obj = sum(y) + 0.0001*sum(sum(sum(sum(r(1:nL,:,:,:)))));
    end
else
    obj = sum(y);
end



ops = sdpsettings('solver','gurobi');
ops.gurobi.TimeLimit = 600;
ops.gurobi.MIPGap = 0.001;

sol = optimize([variablecons, noprovidecons, xx0relationcons,atmostonecons,...
    everyfuncons,x0yrelationcons,nodecapcons,totalcapcons,linearizedcons, ...
    dataratecons, xrzrelationcons, linkcapcons,SFCcons2,...
    mediadelaycons,linkdelaycons,nodedelaycons,E2Econs],obj,ops);
sumy=sum(value(y));
if sumy~=0
    objvalue=value(obj);
else
    objvalue=0;
end

end
