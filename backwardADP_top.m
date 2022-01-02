clearvars
close all
addpath rbfVFA demand_model utilities regtools

% parameters
nc = [5 5 5];   % centers grid 
nw = 1;         % number of workers for parfor loop
nst = 20;       % numner of trajectories for sampling
nss = 20;       % number of state samples at each time
plotFlag = true;
span = 2.0;
threshold = 0; % A threshold on activation, but this is usuallly zero. Let ridge regression take care of ill-conditioning
holdoutPct = 0.25;
widthFactor = 0.3; % \delta in the paper
optModel = "./optimization_model_january"; % optimization model (obtained by running load_optimization_parameters.m)

backwardADP(nc,nw,[nst nss nss],true,span,[holdoutPct threshold],...
    widthFactor,optModel)