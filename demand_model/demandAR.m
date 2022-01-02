function [params,residuals,BIC] = demandAR(traces,na)
% [MODEL,RES,BIC] = DEMANDAR(TRACES,NA) Fits a periodic autoregressive(PAR)
% to the demand trances in TRACES model with constant order NA. The BIC
% is returned for each period and the overall BIC is sum(BIC). See the 
% following references for more on the PAR models:
% 
% [1] McLeod, A. I. (1994). Diagnostic Checking of Periodic Autoregression 
%     Models With Application. Journal of Time Series Analysis, 15(2), 
%     221–233. 
%
% [2] Noakes, D. J., McLeod, A. I., & Hipel, K. W. (1985). Forecasting 
%     monthly riverflow time series. International Journal of Forecasting, 
%     1(2), 179–190. 

tod = traces.tod; % times of day
D = traces.demand; % demand traces
doy = traces.doy; % day of year

Ntod = length(tod);
params.ar = zeros(na+1,Ntod);
params.std = zeros(1,Ntod);
residuals = nan(size(traces.demand));
BIC = nan(1,Ntod);

% initial state distribution from the half hour demand before the first
% tod. Below assumes that the initial state distribution is the same as
% the distribution of midnight states for the whole month.
params.x0.cov = cov( D(end-na+1:end,:)' );
params.x0.mean = mean( D(end-na+1:end,:),2 );
params.tod = tod;

% find the mean of each day of the month
M = length(tod);
mu = mean(D,2,'omitnan');

if na == 0 % mean only model
    params.ar = mu';
    residuals  = D - mu;
    params.std = std(residuals,[],2);
    
    n = length(D);
    np = size(params.ar,1)+length(params.std); % standard deviation estimate too
    BIC = -n/2 * log( params.std.^2 ) + log(n)*np; % Eq. (2.4) in [1]
    return
end

% fit model for first na times of day(requires special construction of 
% regressor matrix X)
DOY = unique(doy);
for ii = 1:na
    X = [];
    Y = [];
    
    % get the previous realization from the same day
    if ii > 1
        b = min([na,ii-1]);
        rowinds1 = ii-b:ii-1;
    else
        b = 0;
        rowinds1 = [];
    end
    b2 = na-b; % realizations from previous day
    rowinds2 = size(D,1)-b2+1:size(D,1);
    rowinds = [rowinds1,rowinds2];
    MUX = mu(rowinds)';
    
    indX = find(doy~=DOY(1)); % exclude first day of month from fitting data
    X = ( D(rowinds,indX-1)'-MUX ); %[TOD,DAY]
    Y = ( D(ii,indX)'-mu(ii) );
    
    p = X\Y;
    params.ar(1:end-1,ii) = p; 
    params.ar(end,ii) = mu(ii) - MUX*p; % bias parameter b (see pooled_stats_and_demand_model.docx)
    r = Y - X*params.ar(1:end-1,ii);
    params.std(ii) = std(r); % find variance for first time
    indFirstDay = (doy==min(doy));
    residuals(ii,~indFirstDay) = r; % residuals for the first na points not possible
    
    n = size(X,1);
    np = length(params.ar(:,ii))+ 1; % +1 for standard deviation estimate
    BIC(ii) = -n/2 * log( params.std(ii)^2 ) + log(n)*np; % Eq. (2.4) in [1]
end

% Fit model to remaining times of day
for kk = na+1:length(tod)
    
    % fit and store AR parameters
    MUX = mu(kk-na:kk-1)';
    X = D(kk-na:kk-1,:)'-MUX; % subtract mean
    Y = (D(kk,:)-mu(kk))';
    p = X\Y;
    params.ar(1:end-1,kk) = p;
    params.ar(end,kk) = mu(kk) - MUX*p; 
    
    % find residual variance
    residuals(kk,:) = Y - X*params.ar(1:end-1,kk);
    params.std(kk) = std(residuals(kk,:));
    
    % calculate BIC for each time period
    n = size(X,1);
    np = length(params.ar(:,kk))+ 1; % include parameter standard deviation
    BIC(kk) = -n/2 * log( params.std(kk)^2 ) + log(n)*np; % Eq. (2.4) in [1]
end

