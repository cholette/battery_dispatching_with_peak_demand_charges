function backwardADP(ncs,nw,Nsx,plotFlag,SPAN,regParam,widthFactor,optModel)

rng('shuffle')
if nw > 1 && isempty(gcp('nocreate'))
    parpool(nw);
end

%% Setup up options
ub = [130;216;130]; % Demand,Ebatt,RecordedPeak
lb = [-60;36;0];
nc = prod(ncs);

%% Set up RBF VFA, samples, and Cost Structure
disp("Width Factor:")
widthFactor
disp("LB:")
lb
disp("UB:")
ub

disp("Loading "+optModel)
dat = load(optModel);
MONTH = month(dat.t0,'shortname');
MONTH = MONTH{1};
sname = "rbfBackwardADP_"+MONTH+"_"+num2str(nc)+".mat";
cons.xub = dat.xub;
cons.xlb = dat.xlb;
cons.ulb = dat.ulb;
cons.uub = dat.uub;
sys = dat.sys;
dd = dat.dd;
T = dat.T;
t0 = dat.t0;
t1 = dat.t1;
t = t0:hours(0.5):t1;
trainModel = dat.trainModel;
disp("Energy cost: "+num2str(sys.energyCost)+", Peak Cost: "+ ...
    num2str(sys.peakCost)+", Damage Cost: "+num2str(sys.damageCost));

w = widthFactor*(ub-lb);
LBRBF = lb - w; % ensure RBF centers are well outside sample domain
UBRBF = ub + w;
sys.outputUpperLimits = inf; % ub(1)+widthFactor*w(1);
sys.outputLowerLimits = -inf; %lb(1)-widthFactor*w(1);
cons.yub = sys.outputUpperLimits;
cons.ylb = sys.outputLowerLimits;

%% State sampling (trajectory)
ns = Nsx(1);
D = zeros(ns,T+1);
rbfVFAFake = createRBFVFA_with_bias(LBRBF,UBRBF,ncs,T,"raw_grid3",SPAN);
sys.epsilon = 1; % all random
parfor jj = 1:ns
    % sample initial state
    sysjj = sys;
    sysjj.x0 = sampleInitialState(sys,trainModel,[cons.xlb;0],[cons.xub;ub(end)],true);
    X = simulatePolicy_RBFVFA(sysjj,dd,rbfVFAFake,t0,t1,[],'q');
    D(jj,:) = X(1,:);
end

[X2,Z]= meshgrid(   linspace(lb(2),ub(2),Nsx(2)),...
                    linspace(lb(3),ub(3),Nsx(3)));
X2 = repmat(X2(:),Nsx(1),1);
Z = repmat(Z(:),Nsx(1),1);
B = ones(Nsx(2)*Nsx(3),1);
Ns = prod(Nsx);

%% Set up RBF VFA, samples, & Cost Structure
rbfVFA = createRBFVFA_with_bias(LBRBF,UBRBF,ncs,T,"raw_grid3",SPAN);
nc = rbfVFA.nc;

%% Backward ADP
ttotal = tic;
for k=T:-1:1 % time loop
    disp("k = "+num2str(k)+"...")
    tloop = tic;
    
    % Get state samples (different for the first time step)
    X1 = kron(D(:,k),B);
    Xsk = [X1';X2'];
    recordedPeakk = Z(:)';
    dk = dd(k);
   
    Vk = zeros(Ns,1);
    disp("Computing value function for "+num2str(Ns)+" samples")
    parfor ii = 1:Ns % Loop over samples
        % Solve Bellman Equation and save optimal value
        xk = Xsk(:,ii);
        recordedPeak = recordedPeakk(ii);
        Vk(ii) = bellmanEquation_RBFVFA(k,xk,recordedPeak,dk,sys,rbfVFA,cons,'quiet');
    end
    disp("Sum of V(:,k): "+num2str(sum(Vk)))
    
    % Estimate RBF coefficients
    X = [Xsk;recordedPeakk];
    [rbfVFA,phis] = estimateRBFParameters(k,X,rbfVFA,Vk,regParam,plotFlag);
    
    if plotFlag
        visualizeTraining(k,t(k),X,phis,Vk,rbfVFA)
    end
    disp(" ... done! Loop time: " + num2str(toc(tloop),'%.2f') + " seconds")
end
ttotal = toc(ttotal);

%% Simulate on test points
x0=[35;38;0]; % demand is mean of jan and aug initial state means
M = 4;
N = 100;
stats = testPolicy_RBFVFA(rbfVFA,trainModel,sys,dd,t0,t1,N,M,true,lb,ub,x0,plotFlag);

% save results
save(sname)
disp("Total time: "+num2str(ttotal))
disp("Model saved as "+sname)
delete(gcp)