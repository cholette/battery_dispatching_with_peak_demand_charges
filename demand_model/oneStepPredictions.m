function Dp = oneStepPredictions(model,traces)

tod = traces.tod; % times of day
D = traces.demand; % demand traces
% doy = traces.doy; % day of year
na = size(model.ar,1)-1;

% compute one-step predictions
Dp = nan(size(D));
for kk = na+1:length(tod)
    
    % fit and store AR parameters
    X = [D(kk-na:kk-1,:);ones(1,size(D,2))];
    Dp(kk,:) = model.ar(:,kk)'*X;
end

