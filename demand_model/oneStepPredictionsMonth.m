function Dp = oneStepPredictionsMonth(model,traces)

t = traces.timestamps; 
D = traces.demand; % demand traces
na = size(model.ar,1)-1;
dt = model.tod(2)-model.tod(1);

% compute one-step predictions
Dp = nan(size(D));
for kk = na+1:length(t)    
    
    todIdx = find( abs(model.tod - timeofday(t(kk))) < 0.1*dt);
    X = [D(kk-na:kk-1,:);ones(1,size(D,2))];
    Dp(kk,:) = model.ar(:,todIdx)'*X;
end