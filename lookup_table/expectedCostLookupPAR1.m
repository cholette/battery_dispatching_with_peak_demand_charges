function [ecost,evalue] = expectedCostLookupPAR1(k,uk,muk,SIGMA,recordedPeak,lookup,sys)

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

% Stage Cost: energy charge (must be in kWh)
alph2 = -muk(end)/sig;
ec = muk(end)*(1-normcdf2(alph2))+normpdf(alph2)*sig;

% Overall stage cost (includes damage)
stageCost = Rd*norm(sys.Dk(:,:,k)*uk) + Re*ec + peakCost*ecdi;

% compute probs of rectangles centered at grid points for the demand and
% state of charge
DgL = xi(:,1)- dx(1)/2;
DgU = xi(:,1)+ dx(1)/2;
SOCgL = xi(:,2)- dx(2)/2;
SOCgU = xi(:,2)+ dx(2)/2;
zgL = xi(:,3)-dx(3)/2;
zgU = xi(:,3)+dx(3)/2;

% adjust boundary points to include infinity
DgU(DgU>max(xi(:,1))) = inf;
SOCgU(SOCgU>max(xi(:,2))) = inf;
zgU(zgU>max(xi(:,3))) = inf;
DgL(DgL<min(xi(:,1))) = -inf;
SOCgL(SOCgL<min(xi(:,2))) = -inf;
zgL(zgL<min(xi(:,3))) = -inf;


SOCkp1 = muk(2);
validSOC = SOCkp1<=SOCgU & SOCkp1>SOCgL;
validPeak = zgU>recordedPeak;
hasPeak = recordedPeak>zgL & recordedPeak<=zgU;
zbar = max(zgL,recordedPeak);
dgL = max(DgL,zbar-sys.Dk(:,:,k)*uk); 
dgU = min(DgU,zgU-sys.Dk(:,:,k)*uk);

m = muk(1);
s = sqrt(SIGMA(1,1));
maskPeakIncrease = validSOC & validPeak & dgU>dgL;
probPeakIncrease = normcdf2(dgU(maskPeakIncrease),m,s) - normcdf2(dgL(maskPeakIncrease),m,s);
maskPeakSame = validSOC & hasPeak & validPeak & dgU>DgL;
probPeakSame = normcdf2(dgU(maskPeakSame),m,s) - normcdf2(DgL(maskPeakSame),m,s);

probs = zeros(1,Ng);
probs(maskPeakIncrease) = probPeakIncrease;
probs(maskPeakSame) = probPeakSame;


evalue = probs*lookup.values(:,k+1);
ecost = stageCost + sys.gamma*evalue;











