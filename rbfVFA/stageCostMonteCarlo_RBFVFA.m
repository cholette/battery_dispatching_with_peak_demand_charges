function [ecost,LB,UB,eStageCost,LBS,UBS] = stageCostMonteCarlo_RBFVFA(k,u,x,recordedPeak,sys,rbfVFA,d,M)

vals = zeros(M,1);
ec = zeros(M,1);
pc = zeros(M,1);
di = zeros(M,1);
dc = sys.damageCost*norm(sys.Dk(:,:,k)*u); % damage cost
for jj = 1:M
    noisew = sqrt(sys.Swk(:,:,k))*randn;
    noisev = sqrt(sys.Svk(:,:,k))*randn;
    yk = sys.Ck(:,:,k)*x+sys.Dk(:,:,k)*u+sys.Ek(:,:,k)*d+sys.Hk(:,:,k)*noisew+noisev;
    xkp1 = [ sys.Ak(:,:,k)*x+sys.Bk(:,:,k)*u+sys.Fk(:,:,k)*d+sys.Gk(:,:,k)*noisew
            yk ];
    xkp1(end) = max([recordedPeak,xkp1(end)]);
    vals(jj) = rbfVFA.theta(:,k+1)'*rbfVFA.phi(xkp1);
    di(jj) = xkp1(end)-recordedPeak;
    
    % stage costs
    ec(jj) = sys.energyCost*(yk>0)*yk; % energy cost
    pc(jj) = sys.peakCost*di(jj); % peak cost
end

stageCosts = dc + pc + ec;
cost =  stageCosts+sys.gamma*vals;
ecost = mean(cost);
eStageCost = mean(stageCosts);


UB = mean(cost)+1.96/sqrt(M)*std(cost);
LB = mean(cost)-1.96/sqrt(M)*std(cost);

UBS = mean(stageCosts)+1.96/sqrt(M)*std(stageCosts);
LBS = mean(stageCosts)-1.96/sqrt(M)*std(stageCosts);