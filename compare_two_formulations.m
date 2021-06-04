
% This file compares the solution efficiency of the two formulations in the
% paper
%
% Optimal Network Slicing for Service-Oriented
% Networks with Flexible Routing and Guaranteed E2E Latency 
%
% by Wei-Kun Chen, Ya-Feng Liu, Antonio De Domenico, Zhi-Quan Luo, and Yu-Hong Dai
% to appear in IEEE Transactions on Network and Service Management, 2021 

% Send your comments and suggestions to   
%             chenweikun@bit.edu.cn      

function comparetwoformulations(iscompactformulation,starttest,endtest)
strname=[num2str(iscompactformulation) '-' num2str(starttest) '-' num2str(endtest)];
objfile=['obj' strname];
soltimefile=['soltime' strname]

fid1=fopen(objfile,'w+');
fid2=fopen(soltimefile,'w+');

for ntest=starttest:endtest
    for k=1:5
        originalgraphfile=['random_networks' '/' 'graph-' num2str(ntest) '.mat'];
        virtualgraphfile=['random_networks' '/' 'virtualgraph' num2str(ntest) '.mat'];
        flowfile=['random_networks' '/' 'SFCs-' num2str(ntest) '.mat'];
        
        if iscompactformulation
            [sol,x,y,r,z,objvalue,sumy,Aver_nodedelay,Aver_linkdelay,...
                SolLinkDelay,SolNodeDelay,averagepath1,maxpath,...
                totalr,averratioonpath,minratioonpath] = ...
                NS_II(2,k,1,1,...
                originalgraphfile,virtualgraphfile,flowfile);
            soltime=sol.solvertime;
        else
            [sol,objvalue] = NS_I(2,k,1,1,...
                originalgraphfile,virtualgraphfile,flowfile);
            soltime=sol.solvertime;
        end
        fprintf(fid1,'%8.4f',objvalue);
        fprintf(fid2,'%12.4f',soltime);
    end
    fprintf(fid1,'\n');
    fprintf(fid2,'\n');
end
fclose(fid1);
fclose(fid2);
end
