% This is the procedure to compute the number of paths and data rates on
% each path in the paper
%
% Optimal Network Slicing for Service-Oriented
% Networks with Flexible Routing and Guaranteed E2E Latency 
%
% by Wei-Kun Chen, Ya-Feng Liu, Antonio De Domenico, Zhi-Quan Luo, and Yu-Hong Dai
% to appear in IEEE Transactions on Network and Service Management, 2021
function [numberofpath, minrateonpath]=minimumpaths(DataRate,nI,nL,L,AllR,originnode,destinctionnode)

maxpath=5;
x=binvar(maxpath,1);
y=sdpvar(maxpath,1);
z=binvar(nL,maxpath,'full');
r=sdpvar(nL,maxpath,'full');

variablecons=[y>=0,r>=0];

%Total data Rate constraint
totaldataratecons=[];
for l=1:nL
    totaldataratecons=totaldataratecons+[sum(r(l,:))==AllR(l)];
end

%Linearized constraint
linearizecons=[];
for l=1:nL
    for p=1:maxpath
        linearizecons=linearizecons+[r(l,p)>=DataRate*z(l,p)+y(p)-DataRate];
        linearizecons=linearizecons+[r(l,p)<=DataRate*z(l,p)];
        linearizecons=linearizecons+[r(l,p)<=y(p)];
    end
end

% x y relation
xycons=[];
for p=1:maxpath
    xycons=xycons+[y(p)<=DataRate*x(p)];
end

% flow conservation constraints
flowconsercons=[];
for i=1:nI
    outgoinglinks=find(L(:,1)==i);
    incominglinks=find(L(:,2)==i);
    if length(outgoinglinks)+length(incominglinks) >0
        addition1=zeros(1,maxpath);
        addition2=zeros(1,maxpath);
        for t=1:length(outgoinglinks)
            addition1=addition1+z(outgoinglinks(t),:);
        end
        
        for t=1:length(incominglinks)
            addition2=addition2+z(incominglinks(t),:);
        end
        if i==originnode
            flowconsercons=flowconsercons+[addition1-addition2==ones(1,maxpath)];
        elseif i==destinctionnode
            flowconsercons=flowconsercons+[addition1-addition2==-1*ones(1,maxpath)];
        else
            flowconsercons=flowconsercons+[addition1-addition2==zeros(1,maxpath)];
        end
    end
end

obj=sum(x);
ops = sdpsettings('solver','gurobi','verbose',0);
ops.gurobi.TimeLimit = 600;
ops.gurobi.MIPGap = 0.001;
sol = optimize([variablecons,totaldataratecons,linearizecons,xycons,flowconsercons],obj,ops);

numberofpath=sum(value(x));

valuey = value(y);
valuex = value(x);
ratios = valuey(valuex>0.999)/DataRate;

minrateonpath = min(ratios);

end

