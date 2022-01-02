clearvars
close all
addpath lookup_table utilities

% parameters
M_ell = 20;
ncs = M_ell*[1 1 1];   	% centers grid 
M_u = 21; 				% control discretization (number of points between -1 and 1)
nw = 6;         		% number of workers for parfor loop
plotFlag = true;
optModel = "./optimization_model_january"; % optimization model (obtained by running load_optimization_parameters.m)
backwardDP_Lookup_PAR1(ncs,nw,plotFlag,M_u,optModel)