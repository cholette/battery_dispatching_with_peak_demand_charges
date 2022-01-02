function x0 = sampleInitialState(sys,demandModel,lb,ub,firstDay,varargin)

na = size(demandModel.ar,1)-1;
x0 = zeros(na+2,1);
x0(1:na) = mvnrnd(demandModel.x0.mean,demandModel.x0.cov); % Initial Demand
if nargin > 5 % set up beta distribution parameters for recorded peak and SOC sampling
    aSOC = varargin{1}(1,1);
    bSOC = varargin{1}(1,2);
    aPeak = varargin{1}(2,1);
    bPeak = varargin{1}(2,2);
else
    aSOC = 1;
    bSOC = 1;
    aPeak = 1;
    bPeak = 1;
end

x0(na+1) = sys.stateLowerLimits(2)+...
    (sys.stateUpperLimits(2)-sys.stateLowerLimits(2))*betarnd(aSOC,bSOC); %SoC

if firstDay
    recPeak = 0; % recorded peak is zero on first day
else
    %             x0(3) = lb(3) + (ub(3)-lb(3))*rand;
    recPeak = lb(3) + (ub(3)-lb(3))*betarnd(aPeak,bPeak);
end
x0(na+2) = recPeak;