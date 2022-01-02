% Script for setting up the state space model and other parameters needed
% for the backward ADP algorithm.

addpath utilities
clearvars, close all
mm = 1; % Number of the month of interest. 1 and 8 used in the paper
yyyy = 2017;
DAYS = 31;
peakCost = 22.436; % 22.436; %$/kW peak
damageCost = 0.01; % $/kW battery damage cost of charge/discharge;
energyCost = 14.961/2/100; % $/kW/half hour. From 14.961 c/kWh 
gam = 1; % discount rate. Please note that only gam=1 has been tested
months = ["january","february","march","april","may","june","july","august",...
    "september","october","november","december"];
dfile = "./demandModel_"+months(mm); % demand model (obtained from fit_demand_model.m)
name = "./optimization_model_"+months(mm); % save name

%% setup optimization problem
T = 48*DAYS;
t0 = datetime(yyyy,mm,1); % beginning of month
dat = load(dfile);
traces = dat.tracesTrain;
t1 = t0 + hours(T*0.5-0.5);
t = t0:hours(0.5):t1;

% base system (nans are placeholders for the time-varying parameters)
base_system.A = [nan,0;0,1];
base_system.B = [0;18.6];
base_system.C = [nan,0];
base_system.D = [40];
base_system.F = [1;0];
base_system.E = 1;
base_system.G = [1;0];
base_system.H = 1;
[sys,dd] = buildSystemMatrices(base_system,dat.trainModel,t); % use the training model (which excludes some data and is only AR(1))
sys.dt = duration(0.5,0,0); % time step in hours

% cost function parameters
sys.energyCost = energyCost;
sys.damageCost = damageCost;
sys.peakCost = peakCost;
sys.gamma = gam; % gam = 1 for the paper. gam < 1 untested. 

% constraints
sys.stateLowerLimits = [-inf;36];
sys.stateUpperLimits = [inf;216];
sys.controlUpperLimits = [1]; % -1<=u<=1, % of max power (max power absorbed into B)
sys.controlLowerLimits = [-1];
xub = sys.stateUpperLimits;
xlb = sys.stateLowerLimits;
uub = sys.controlUpperLimits;
ulb = sys.controlLowerLimits;

% training traces and model
trainModel = dat.trainModel;
tracesTrain = dat.tracesTrain;

% test traces and model
testModel = dat.testModel;
tracesTest = dat.tracesTest;

save(name)