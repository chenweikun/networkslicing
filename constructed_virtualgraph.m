
% This file is to demonstrate how to transform the physical network into the
% virtual networks  in the paper 
%
% Optimal Network Slicing for Service-Oriented
% Networks with Flexible Routing and Guaranteed E2E Latency 
%
% by Wei-Kun Chen, Ya-Feng Liu, Antonio De Domenico, Zhi-Quan Luo, and Yu-Hong Dai
% to appear in IEEE Transactions on Network and Service Management, 2021 

% Send your comments and suggestions to   
%             chenweikun@bit.edu.cn        


function constructed_virtualgraph(nI,nS,nV,nL,V,L,Link_cap,Link_delay,filename)

nnewI=nI+(nS)*nV;
nnewV=nV*nS;
nnewL=(nS+1)*nV + nL;
newI = 1:nnewI;
newV = nI+1:nnewI;

newL=L;
newLink_cap=Link_cap;
newLink_delay=Link_delay;


k=nL+1;
for i=1:nV
    for j=(nI+(i-1)*nS+1):(nI+i*nS)
        newL(k,1)=V(i);
        newL(k,2)=j;
        newLink_cap(k)=100000;
        newLink_delay(k)=0;
        k=k+1;
        newL(k,1)=j;
        newL(k,2)=V(i);
        newLink_cap(k)=100000;
        newLink_delay(k)=0;
        k=k+1;
    end
end

save(filename,'newL','newLink_cap','newLink_delay');


end
