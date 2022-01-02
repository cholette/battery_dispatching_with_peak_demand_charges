function [ecost,evalue] = expectedCostLookup(k,uk,muk,SIGMA,recordedPeak,lookup,sys) %#codegen

Ng = lookup.Ng;
Re = sys.energyCost;
Rd = sys.damageCost;
peakCost = sys.peakCost;
xi = lookup.xi; % last one is peak
dx = lookup.dx;
na = size(lookup.xi,2)-2;

muk = muk + [sys.Bk(:,:,k);sys.Dk(:,:,k)]*uk;

if recordedPeak < 0
    recordedPeak = 0;
end

sig = sqrt(SIGMA(end,end));
D = muk(end)-recordedPeak;
alph = -D/sig;

% Stage cost: peak increase
ecdi = D*(1-normcdf2(alph))+sig*normpdf(alph);

% Stage Cost: energy charge (cost in $/kWh)
alph2 = -muk(end)/sig;
ec = muk(end)*(1-normcdf2(alph2))+normpdf(alph2)*sig;

% Overall stage cost (includes damage)
stageCost = Rd*norm(sys.Dk(:,:,k)*uk) + Re*ec + peakCost*ecdi;

% compute probs of rectangles centered at grid points for the demand and
% state of charge
lb = xi-dx/2;
ub = xi+dx/2;
includesPeak = (recordedPeak<ub(:,end) & recordedPeak>=lb(:,end) );
includesSOC = ( muk(na+1)<ub(:,na+1) & muk(na+1)>=lb(:,na+1) );
validPeak = (ub(:,end)>=recordedPeak);
validDemand = (lb(:,1)<=ub(:,end));
includesAll = includesSOC & includesPeak & validPeak & validDemand;

xg = xi(:,1:end-1);
zg = xi(:,end);
zk = recordedPeak*ones(Ng,1);
zlb1 = max([zg-dx(end)/2,zk],[],2);
zlb2 = -inf(sum(includesAll),1);
zub1 = max([zg+dx(end)/2,zlb1+1e-6],[],2);
zub2 = min([zg+dx(end)/2,zk],[],2);
xlb = xg-dx(1:end-1)/2;
xub = xg+dx(1:end-1)/2;

xlb(xlb<min(xg)) = -inf;
xub(xub>max(xg)) = inf;
zlb1(zlb1<min(zg)) = -inf;
zlb2(zlb2<max(zg)) = -inf;
zub1(zub1>max(zg)) = inf;
zub2(zub2>max(zg)) = inf;

mask1 = (includesSOC & validPeak & validDemand);
I1t = mvncdf([xlb(mask1,:),zlb1(mask1,:)],...
    [xub(mask1,:),zub1(mask1,:)],muk',SIGMA);
I1 = zeros(Ng,1);
I1(mask1) = I1t;

I2t = mvncdf([xlb(includesAll,:),zlb2],[xub(includesAll,:),...
    zub2(includesAll)],muk',SIGMA);
I2 = zeros(Ng,1);
I2(includesAll) = I2t;

% deal with NaNs
I1(isnan(I1))=0;
I2(isnan(I2))=0;
probs = I1 + I2;


evalue = probs'*lookup.values(:,k+1);
ecost = stageCost + sys.gamma*evalue;











