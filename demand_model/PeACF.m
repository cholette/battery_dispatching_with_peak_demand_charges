function [lags,rho,mu,bnds] = PeACF(Y,Nlags,S)
% Perodic autorrelation function (can have multiple realizations of a time
% series) with perodicity S. Rows are time samples and columns are different
% realizations. See [1] for more.
%
% [1]   A. I. McLeod, “Diagnostic Checking of Periodic Autoregression Models With Application,”
%       J. Time Ser. Anal., vol. 15, no. 2, pp. 221–233, 1994.

% build lags vector
lags = 1:Nlags;
N = size(Y,1); % length of each series
M = size(Y,2); % number of time series
rho = zeros(Nlags+1,S); % autocovariance matrix
na = nan(Nlags+1,S); % number of samples available (should be n for all but first few m)
n = N*M/S; % number of repetitions of a day over all time series
Z = [];
for ii = 1:M
    Z = [Z,reshape(Y(:,ii),S,[])];
end
% Z = reshape(Y',n,S); % reshape to have all samples with same periodic time index in same column
mu = mean(Z,2,'omitnan')'; % mean of each sample
c0 = var(Z,1,2,'omitnan');
rho(1,:) = 1;
for m=1:S % periodic index
    for k = lags
        ind = m:S:size(Y,1); % vector of samples with same periodic index
        indLag = ind-k; % vector of samples with same lag from periodic index m
        
        ind(indLag<=0)=[]; % remove the negative indices
        indLag(indLag<=0)=[];
        na(k+1,m) = length(indLag)*M;
        
        D = (Y(ind,:)-mu(m))/sqrt(c0(m)); % normalize so I get an Autocorrelation instead of Autocovariance
        if m-k>0
            indmu = m-k;
        elseif m-k == 0
            indmu = m;
        else
            indmu = S+(m-k)+1;
        end
        DL = (Y(indLag,:)-mu(indmu))/sqrt(c0(indmu));
        
        rho(k+1,m) = 1/na(k+1,m)*sum(sum(D.*DL,'omitnan'));
    end
end

lags = [0,lags];
bnds = [-1.96/sqrt(n) 1.96/sqrt(n)];