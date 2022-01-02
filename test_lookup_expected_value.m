% % test expected cost for lookup tables
clearvars, close all
addpath lookup_table utilities
rbfVFA_file = "./LookupBackwardDP_21_8000_1";
dat = load(rbfVFA_file); % lookup table VFA

% state
N = size(dat.sys.Ak,3);
M = 10000;
tol = 1e-6;
k = randi(N);
xk =  [100*rand;50+150*rand]; %[150;216];
uk = -1+2*rand;
dk = dat.dd(k); 
% dk = dat.dd(k); 
recordedPeak = max([100*rand,xk(1)]);
sys = dat.sys;
sys.gamma = 1;

% formula
% lookupVFA = dat.lookupTable;
lookupVFA = createLookupTable(20*ones(1,3),[-50,0,0],[150,216,200],N);
lookupVFA.values = 1000*rand(lookupVFA.Ng,N);
% [ecost,evalue] = expectedCostLookup(k,uk,xk,recordedPeak,dk,lookupVFA,sys);
[muk,S] = extractSystemInfo(sys,xk,dk,k,tol);

t2 = tic;
% [ecost2,evalue2] = expectedCostLookup(k,uk,muk,S,recordedPeak,lookupVFA,sys,M);
[ecost2,evalue2] = expectedCostLookup(k,uk,muk,S,recordedPeak,lookupVFA,sys);
runTime2 = toc(t2);

t3 = tic;
[ecost3,evalue3] = expectedCostLookupPAR1(k,uk,muk,S,recordedPeak,lookupVFA,sys);
runTime3 = toc(t3);

% monte carlo
t4 = tic;
w = mvnrnd(0,sys.Swk(:,:,k),M);
v = mvnrnd(0,sys.Svk(:,:,k),M);
COSTS = nan(1,M);
COSTT = nan(1,M);
BMU = nan(1,M);
for mm = 1:M
    xkp1 = sys.Ak(:,:,k)*xk+sys.Bk(:,:,k)*uk+sys.Fk(:,:,k)*dk + sys.Gk(:,:,k)*w(mm);
    yk = sys.Ck(:,:,k)*xk+sys.Dk(:,:,k)*uk+sys.Ek(:,:,k)*dk+sys.Hk(:,:,k)*w(mm)+v(mm);
    zk = max([yk,recordedPeak]);
    stageCost = sys.damageCost*norm(sys.Dk(:,:,k)*uk) + sys.energyCost*(yk>0)*yk +...
        sys.peakCost*(zk-recordedPeak);
    [~,indg] = aggragationFunction([xkp1;zk]',lookupVFA,'edges');
    BMU(mm) = indg;
    COSTS(mm) = stageCost;
    COSTT(mm) = sys.gamma*lookupVFA.values(indg,k+1);
end
COST = COSTS+COSTT;
sterr = std(COST)/sqrt(M);
mu = mean(COST);
runTime4 = toc(t4);

disp("Mean Cost 95% CI (Sim): ["+num2str(mu-1.96*sterr,'%.2f')+", "+num2str(mu+1.96*sterr,'%.2f')+"]"+...
    " in "+num2str(runTime4*1000)+" ms")
disp("Expected Cost (MVNCDF): "+num2str(ecost2,'%.2f')+" in "+num2str(runTime2*1000)+" ms");
disp("Expected Cost (NORMCDF2): "+num2str(ecost3,'%.2f')+" in "+num2str(runTime3*1000)+" ms")