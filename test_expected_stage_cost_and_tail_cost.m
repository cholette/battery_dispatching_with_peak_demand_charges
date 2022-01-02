% The below script was used to verify the expected value computation for the 
% Bellman equation when employing an RBF VFA (Theorem 2) with Monte 
% Carlo simulations. 

clearvars, close all
addpath utilities rbfVFA

% Example models
load ./rbfBackwardADP_Jan_125_1 % RBF VFA .mat file (obtaind by running backwardADP_top.m)
load ./demandModel_January % demand model (obtained from fit_demand_model.m)

toStr = @(x) num2str(x,'%.4f');
inRange = @(x,lb,ub) x<=ub & x>=lb;
tol = 1e-4; % tolerance condition number for covariance matrix (\varepsilon in paper)
M = 1e4; % number of Monte Carlo simulations for expected value estimation
numSims = 100; % number of different initial conditions used to test computations

% set up verification results structure
results.stageCosts.test = nan(numSims,1);
results.stageCosts.delta = nan(numSims,1);
results.stageCosts.simulated = nan(numSims,1);
results.tailCosts.test = nan(numSims,1);
results.tailCosts.delta = nan(numSims,1);
results.tailCosts.simulated = nan(numSims,1);
results.totalCosts.test = nan(numSims,1);
results.totalCosts.delta = nan(numSims,1);
results.totalCosts.simulated = nan(numSims,1);

% cost structure (random peak, energy and damage costs)
sys.peakCost = 100*rand;
sys.energyCost = 1000*rand;
sys.damageCost = 10*rand;

for nn = 1:numSims
    ii = randi(size(sys.Swk,3)); % sample time index
    u = sys.controlLowerLimits+rand*(sys.controlUpperLimits-...
        sys.controlLowerLimits); % sample input
    dk = -10+100*rand; % sample disturbance
    xk = lb(1:2)+rand(2,1).*(ub(1:2)-lb(1:2)); % sample state
    ubr = max([ub(end),xk(1)]);
    recordedPeak = lb(end)+rand*(ubr-lb(end)); % sample peak;
    
    disp("ii="+num2str(ii)+", xk = ["+toStr(xk(1))+","+toStr(xk(2))+...
        "]; recordedPeak = "+toStr(recordedPeak))
    
    muk = [ sys.Ak(:,:,ii)*xk+sys.Bk(:,:,ii)*u+sys.Fk(:,:,ii)*dk
        sys.Ck(:,:,ii)*xk+sys.Dk(:,:,ii)*u+sys.Ek(:,:,ii)*dk];
    SIGMA = [   sys.Gk(:,:,ii)*sys.Swk(:,:,ii)*sys.Gk(:,:,ii)'      sys.Gk(:,:,ii)*sys.Swk(:,:,ii)*sys.Hk(:,:,ii)'
        sys.Hk(:,:,ii)*sys.Swk(:,:,ii)*sys.Gk(:,:,ii)'      sys.Hk(:,:,ii)*sys.Swk(:,:,ii)*sys.Hk(:,:,ii)' + sys.Svk(:,:,ii) ];
    SIGMAs = SIGMA;
    
    
    %% monte carlo (vals is mvn while vals2 is one simulation at a time)
    tic
    xp = mvnrnd(muk',SIGMAs,M);
    xp(:,end) = max([xp(:,end),recordedPeak*ones(M,1)],[],2);
    vals = zeros(M,1);
    vals2 = zeros(M,1);
    for jj = 1:M
        vals(jj) = rbfVFA.theta(:,ii+1)'*rbfVFA.phi(xp(jj,:)');
        noisew = sqrt(sys.Swk(:,:,ii))*randn;
        noisev = sqrt(sys.Svk(:,:,ii))*randn;
        xkp1 = [ sys.Ak(:,:,ii)*xk+sys.Bk(:,:,ii)*u+sys.Fk(:,:,ii)*dk+sys.Gk(:,:,ii)*noisew
            sys.Ck(:,:,ii)*xk+sys.Dk(:,:,ii)*u+sys.Ek(:,:,ii)*dk+sys.Hk(:,:,ii)*noisew+noisev];
        xkp1(end) = max([recordedPeak,xkp1(end)]);
        vals2(jj) = rbfVFA.theta(:,ii+1)'*rbfVFA.phi(xkp1);
    end
    EVkp1Sim = mean(vals2); 
    UCB = EVkp1Sim +1.96/sqrt(M)*std(vals2);
    LCB = EVkp1Sim -1.96/sqrt(M)*std(vals2);   
    
    %% Check expected cost function (uses method from Theorem 2)
    [mukTest,SIGMAkTest] = extractSystemInfo(sys,xk,dk,ii,tol);
    tic
    [ecost,EVkp1,stageCost] = expectedCost_RBFVFA(ii,u,xk,recordedPeak,sys,...
        rbfVFA,mukTest,SIGMAkTest);
    
    % test tail costs
    disp("")
    disp("=========== Monte Carlo tail cost vs Analytical tail cost ==============")
    disp("Monte Carlo: ["+toStr(LCB) + " " + toStr(UCB)+"]"+" in "+num2str(toc)+" sec")
    disp("Analytical: "+toStr(EVkp1)+" in "+num2str(toc)+" sec")
    results.tailCosts.simulated(nn) = EVkp1Sim;
    results.tailCosts.delta(nn) = EVkp1Sim - EVkp1;
    if EVkp1 <= UCB+eps && EVkp1 >= LCB-eps
        disp('...Pass!')
        results.tailCosts.test(nn) = 1;
    elseif EVkp1 < UCB-eps
        disp("Fail. stageCost is "+num2str(LCB-EVkp1,'%.2e')+" below lower CI")
        results.tailCost.test(nn) = 0;
    elseif EVkp1 > UCB + eps
        disp("Fail. stageCost is "+num2str(EVkp1-UCB,'%.2e')+" above upper CI")
        results.tailCost.test(nn) = 0;
    elseif EVkp1 > UCB + eps
        disp("Fail. stageCost is "+num2str(EVkp1-UBS,'%.2e')+" above upper CI")
        results.tailCost.test(nn) = 0;
    end
    
    [ecostMC,LBMC,UBMC,stageCostMC,LBS,UBS] = stageCostMonteCarlo_RBFVFA(ii,u,xk,...
        recordedPeak,sys,rbfVFA,dk,M);
    
    disp("")
    disp("")
    disp("================= Monte Carlo stage cost vs Analytical stage cost ======================")
    disp("Monte Carlo: ["+toStr(LBS)+","+toStr(stageCostMC)+","+toStr(UBS)+"]")
    disp("Analytical: "+toStr(stageCost))
    results.stageCosts.simulated = stageCostMC;
    results.stageCosts.delta(nn) = stageCostMC - stageCost;
    if stageCost <= UBS+eps && stageCost >= LBS-eps
        disp('...Pass!')
        results.stageCosts.test(nn) = 1;
    elseif stageCost < LBS-eps
        disp("Fail. stageCost is "+num2str(LBS-stageCost,'%.2e')+" below lower CI")
        results.stageCosts.test(nn) = 0;
    elseif stageCost > UBS + eps
        disp("Fail. stageCost is "+num2str(stageCost-UBS,'%.2e')+" above upper CI")
        results.stageCosts.test(nn) = 0;
    end
    
    disp("")
    disp("================= Monte Carlo total cost vs Analytical total cost ==================")
    disp("Monte Carlo: [" + toStr(LBMC)+", "+toStr(ecostMC)+", "+toStr(UBMC)+']')
    disp("Analytical: " + toStr(ecost))
    results.totalCosts.simulated(nn) = ecostMC;
    results.totalCosts.delta(nn) = ecostMC - ecost;
    if ecost <= UBMC+eps && ecost >= LBMC-eps
        disp('...Pass!')
        results.totalCosts.test(nn) = 1;
    elseif ecost < LBMC-eps
        disp("Fail. Total cost is "+num2str(LBMC-ecost,'%.2e')+" below lower CI")
        results.totalCosts.test(nn) = 0;
    elseif ecost > LBMC + eps
        disp("Fail. Total cost is "+num2str(ecost-UBMC,'%.2e')+" above upper CI")
        results.totalCosts.test(nn) = 0;
    end
end

%% Plots
figure
histogram(results.totalCosts.delta./results.totalCosts.simulated)
xlabel("relative error (%)")
ylabel("Probability")
title("Relative error of simulation vs. analytical formula (N="+...
    num2str(numSims,'%d')+")")
h = gca;
h.FontSize = 14;

figure
histogram(results.totalCosts.test,'normalization','probability')
h = gca;
h.FontSize = 14;
h.XTick = [0 1];
h.XTickLabel = ["Outside CI","Within CI"];
h.XTickLabelRotation = 0;
hold on
plot(h.XLim,0.05*[1 1],'r--','linewidth',2)
ylabel("Probability")
title("% of analytical total costs within"+newline+"95% CI of simulation mean")
text(0,0.05,"5% limit",'Fontsize',12,'Color','r','HorizontalAlignment',...
    'center','VerticalAlignment','bottom')