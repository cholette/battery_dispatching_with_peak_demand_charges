function stats = testPolicy_RBFVFA(rbfVFA,trainModel,sys,dd,t0,t1,N,M,firstDay,lb,ub,varargin)

if nargin >= 12
    x0Test = varargin{1};
    plotFlag = varargin{2};
else
    x0Test = [];
    plotFlag = true;
end


mu = trainModel.x0.mean;
sig = sqrt(trainModel.x0.cov);
if firstDay
    M2 = floor(M^(1/2));
    dGrid = linspace(mu-2*sig,mu+2*sig,M2);
    socGrid = linspace(lb(2),ub(2),M2);
    [dG,socG] = meshgrid(dGrid,socGrid);
    peakG = zeros(size(dG));
    X0 = [dG(:)';socG(:)';peakG(:)'];
else
    M2 = floor(M^(1/3));
    dGrid = linspace(mu-3*sig,mu+3*sig,M2);
    socGrid = linspace(lb(2),ub(2),M2);
    peakGrid = linspace(lb(3),ub(3),M2);
    [dG,socG,peakG] = meshgrid(dGrid,socGrid,peakGrid);
    X0 = [dG(:)';socG(:)';peakG(:)'];
end

% add the test point if it exists
if ~isempty(x0Test)
    X0 = [X0,x0Test];
end

% simulate
sys.epsilon = 0; % greedy
stats = evaluatePolicy_RBFVFA(rbfVFA,sys,dd,t0,t1,X0,N);

if plotFlag % plot results of test points
    figure
    boxchart(stats.cost)
    hold on
    plot(stats.estCost,'g*')
    xlabel('Inital State index')
    ylabel('Cost')
    
    figure
    scatter3(X0(1,:),X0(2,:),X0(3,:),[],mean(stats.cost,1))
    xlabel('Demand')
    ylabel('SOC')
    zlabel('Recorded Peak')
    
    figure
    X = repmat(stats.estCost,size(stats.cost,1),1); 
    Y = stats.cost; 
    scatter(X(:),Y(:));
    hold on
    h = gca;
    xlim([0,h.XLim(2)]);
    plot([0,h.XLim(2)],[0,h.XLim(2)],'r--') % ideal line
    xlabel('RBF prediction')
    ylabel('Cost')
    
    
end

% print results for a test point
if isempty(x0Test) % pick a point
    ind = find(X0(1,:)==dGrid(ceil(M2/2)) & X0(2,:)==socGrid(1) & ...
        X0(3,:)==peakGrid(1));   
else
    ind = size(X0,2);
end

COST = stats.cost(:,ind);
PRED = stats.estCost(ind);
disp('================ Summary ===========================')
UB = mean(COST)+1.96/sqrt(N)*std(COST);
LB = mean(COST)-1.96/sqrt(N)*std(COST);
disp(['Test Initial State: [',num2str(X0(:,ind)','%.2f,%.2f,%.2f'),']'])
disp(['Controlled Cost: [',num2str(LB,'%.2f'),', ',num2str(UB,'%.2f'),']'])
disp(['Predicted Cost (RBF): ',num2str(PRED,'%.2f')])