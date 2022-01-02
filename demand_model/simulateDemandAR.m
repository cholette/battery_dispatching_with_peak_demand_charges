function [D,Xp] = simulateDemandAR(params,times,varargin)
% times must be a datetimes grid

N = length(times);
na = size(params.ar,1)-1;
tod = params.tod;
D = zeros(N,na);
Xp = zeros(N,na);

dt = diff(times);
if ~all( abs(dt-dt(1))<1e-10 )
    error('Time vector needs to be evenly sampled')
else
    dt = dt(1);
end

% generate from initial state distribution or start from known x0
if nargin <=2
    D(1,1:na) = mvnrnd(params.x0.mean,params.x0.cov); % sample initial distribution
%     t(1) = datetime(traces.year(1),1,1) + days(traces.doy(1)-1) + traces.tod(1); % start at beginning of traces
    Xp(1,1:na) = params.x0.mean;
else
    D(1,1:na) = varargin{1};
    Xp(1,1:na) = D(1:na);
end

% simulate
for kk = 2:N
    if class(times(kk-1))~="duration"
        nowtod = timeofday(times(kk-1));
    else
        nowtod = times(kk-1);
    end
%     t(kk) = t(kk-1)+traces.dt;
    ind = find(tod==nowtod);
    
    % forecast and simulate
    Xp(kk,1) = [Xp(kk-1,:),1]*params.ar(:,ind);
    Xp(kk,2:end) = Xp(kk-1,1:end-1);
    D(kk,1) = [D(kk-1,:),1]*params.ar(:,ind) + params.std(ind)*randn;
    D(kk,2:end) = D(kk-1,1:end-1);
end



