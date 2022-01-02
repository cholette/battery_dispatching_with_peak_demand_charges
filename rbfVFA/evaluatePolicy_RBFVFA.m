function stats = evaluatePolicy_RBFVFA(rbfVFA,sys,dd,t0,t1,X0,N,varargin)

if nargin<7 || isempty(varargin)
    verbose = true;
else
    verbose = varargin{1};
end

M = size(X0,2); % number of evaluation points
costUC = nan(N,M);
cost = nan(N,M);
estCost = nan(1,M);
controlledPeak = nan(N,M);
uncontrolledPeak = nan(N,M);
flag = iscell(rbfVFA.phi);
parfor mm = 1:M
    sysmm = sys;
    sysmm.x0 = X0(:,mm);
    
    if flag
        estCost(mm) = rbfVFA.theta{1}'*rbfVFA.phi{1}(sysmm.x0);
    else
        estCost(mm) = rbfVFA.theta(:,1)'*rbfVFA.phi(sysmm.x0);
    end
    
    for nn = 1:N
        
        % simulate
        [X1,~,COST,Y] = simulatePolicy_RBFVFA(sysmm,dd,rbfVFA,t0,t1,[],"quiet");
        
        % get stats
        cost(nn,mm) = sum(COST);
        [cdMax,~] = max(Y);
        [uncontrolledPeak(nn,mm),~] = max(X1(1,2:end));
        controlledPeak(nn,mm) = max([sysmm.x0(end),cdMax]);
        
        YUC = X1(1,2:end,1);
        mYUC = max([sysmm.x0(end),YUC]);
        costUC(nn,mm) = dot(sysmm.energyCost*(YUC>0),YUC) + ...
        mYUC*sysmm.peakCost; % only works if gamma is one!
        
        if verbose
            disp(['x0= [',num2str(sysmm.x0','%.2f,%.2f,%.2f'),'], m=',num2str(mm),', n= ',num2str(nn),', Controlled (Sim): ',num2str(cost(nn,mm),'%.0f'),...
                ', Predicted Expected Cost (RBF): ',num2str(estCost(mm),'%.2f')])
        end
        
    end
end

% pack up stats
stats.costUncontrolled = costUC;
stats.cost = cost;
stats.estCost = estCost;
stats.controlledPeak = controlledPeak;
stats.uncontrolledPeak = uncontrolledPeak;
stats.x0 = X0;