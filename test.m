

% This contains a simple example for the network slicing problem in the paper
%
% Optimal Network Slicing for Service-Oriented
% Networks with Flexible Routing and Guaranteed E2E Latency 
%
% by Wei-Kun Chen, Ya-Feng Liu, Antonio De Domenico, Zhi-Quan Luo, and Yu-Hong Dai
% to appear in IEEE Transactions on Network and Service Management, 2021

% Send your comments and suggestions to   
%             chenweikun@bit.edu.cn        

nP=2;
nK1=2;
ntest=1;
isdelaycons=1;
isdelayobj=1;


originalgraphfile='random_networks/exgraph.mat';
virtualgraphfile='random_networks/exvirtualgraph.mat';
flowfile='random_networks/exSFCs.mat';

% run the compact (NS-II) formulation
[sol,x,y,r,z,objvalue,sumy,Aver_nodedelay,Aver_linkdelay,...
    SolLinkDelay,SolNodeDelay,averagepath,maxpath] = ...
    NS_II(nP,nK1,isdelaycons,isdelayobj,...
    originalgraphfile,virtualgraphfile,flowfile);

% run the (NS-I) formulation
[sol,objvalue] = NS_I(nP,nK1,isdelaycons,isdelayobj,...
    originalgraphfile,virtualgraphfile,flowfile);


