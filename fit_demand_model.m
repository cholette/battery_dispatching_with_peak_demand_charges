clearvars
addpath demand_model utilities

% Script for identifying demand model. Note that leap year handling won't
% work for February (and neither will the below script).
% Also, while NA>1 will work in this script, we will not be able to convert
% it to state space using the NA > 1 the state space function conversion
% in load_system_matrices.m won't work because it was written only for 
% the NA=1 special case (which was the case used in the paper)

naTrain = 1; % Autoregressive order for training model
naTest = naTrain; % Autoregressive order for test model
maxLags = 1*48; % Maximum lags for PACF plotting
mm = 1; % month
includeTrainingInTest = true;
percentSigmaToAdd = 0;
yyyy = 2016; % year for end of training data. Data was present from July 2013 to April 2017 in the paper.
dataFile = 'data/demand_data.xlsx'; % This data file is not provided due to IP arrangements. See the README.md for the format of this file
dataSheet = 'data'; % Sheet name where the demand data is stored

t0 = datetime(yyyy,mm,1,0,0,0); % beginning of m
EOM = eomday(yyyy,month(t0));
endTrain = datetime(yyyy,mm+1,1,0,0,0,0); % end of training data
if leapyear(yyyy) && mm > 2 % if it's a leap year, need to subtract 1 of doy after Feb
    doy1 = day(t0,'dayofyear')-1;
    doy2 = doy1+ EOM-1;
else
    doy1 = day(t0,'dayofyear');
    doy2 = doy1+ EOM-1;
end

mm = month(t0,'name');
savename = ['./demandModel_',mm{1}];

% load in data, split into train and test demand profiles
[datTrain,datTest] = loadData(dataFile,dataSheet,endTrain,includeTrainingInTest);
tracesTrain = demandProfiles(datTrain,doy1,doy2,'group');
tracesTest = demandProfiles(datTest,doy1,doy2,'group');

% fit PAR models
[trainModel,trainResiduals,trainBIC] = demandAR(tracesTrain,naTrain);
[testModel,testResiduals,testBIC] = demandAR(tracesTest,naTest);

% add variance to training model
trainModel.std = (1+percentSigmaToAdd)*trainModel.std;
save(savename);

%% Analyse residuals on TRAINING data
figure("Name","Training Model")
Ny = length(tracesTrain.year);
% [lags,r,bnds] = autocorr2(reshape(trainResiduals,[EOM*48,Ny]),48*5); 
tracesTrainMonth = demandProfilesMonth(datTrain,doy1,doy2);
Dp = oneStepPredictionsMonth(trainModel,tracesTrainMonth);
trainResiduals2 = tracesTrainMonth.demand - Dp;
[lags,r,mu,bnds] = PeACF(trainResiduals2,maxLags,48);
stem(lags,r); 
h = gca; hold on; 
plot(h.XLim,[bnds(1) bnds(1)],'g--');
plot(h.XLim,[bnds(2) bnds(2)],'g--');
xlabel('lags')
ylabel('PeACF')

%% Analyse residuals on TESTING data
figure("Name","Testing Model")
Ny = length(tracesTest.year);
% [lags,r,bnds] = autocorr2(reshape(testResiduals,[EOM*48,Ny]),48*5); 
tracesTestMonth = demandProfilesMonth(datTest,doy1,doy2);
Dp = oneStepPredictionsMonth(testModel,tracesTestMonth);
testResiduals2 = tracesTestMonth.demand - Dp;
[lags,r,mu,bnds] = PeACF(testResiduals2,maxLags,48);
stem(lags,r); 
h = gca; hold on; 
plot(h.XLim,[bnds(1) bnds(1)],'g--');
plot(h.XLim,[bnds(2) bnds(2)],'g--');
xlabel('lags')
ylabel('PeACF')








