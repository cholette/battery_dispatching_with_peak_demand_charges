function [ecost,EVkp1,stageCost] = expectedCost_RBFVFA(k,uk,xk,recordedPeak,sys,rbfVFA,muk,SIGMA)

% unpack coefficients
if iscell(rbfVFA.xi)
    ncoeff = size(rbfVFA.theta{k+1},1);
    nc = rbfVFA.nc(k+1);
    if ncoeff == nc
        biasFlag = false;
    elseif ncoeff == nc+1
        biasFlag = true;
    else
        error('Number of coefficients in the VFA is inconsistent with the number of centers')
    end
    
    cm = rbfVFA.xi{k+1};
    sigwidth = rbfVFA.sigma(:,:,k+1);
    wm = rbfVFA.theta{k+1}(1:nc);
    if biasFlag
        wmb = rbfVFA.theta{k+1}(nc+1);
    else
        wmb = 0;
    end
else
    ncoeff = size(rbfVFA.theta,1);
    nc = rbfVFA.nc;
    if ncoeff == nc
        biasFlag = false;
    elseif ncoeff == nc+1
        biasFlag = true;
    else
        error('Number of coefficients in the VFA is inconsistent with the number of centers')
    end
    cm = rbfVFA.xi;
    sigwidth = rbfVFA.sigma;
    wm = rbfVFA.theta(1:nc,k+1);
    if biasFlag
        wmb = rbfVFA.theta(nc+1,k+1);
    else
        wmb = 0;
    end
end

% handle variable sigma
% sig = nan;
if length(sigwidth) > 1
    %     sig = rbfVFA.sigma(k+1);
%     GAMMA = diag(sigwidth.^2);
    GAMMAi = diag(1./(sigwidth.^2));
else
    %     sig = rbfVFA.sigma;
%     GAMMA = sigwidth*eye(length(xk)+1);
    GAMMAi = 1./(sigwidth.^2)*eye(length(xk)+1);
end


Re = sys.energyCost;
Rd = sys.damageCost;
peakCost = sys.peakCost;

muk = muk + [sys.Bk(:,:,k)*uk;sys.Dk(:,:,k)*uk];

% Stage cost: recorded peak component
% sig = sqrt(Hk*Swk*Hk'+Svk);
% mu = Ck*xk + Ek*dk + Dk*uk; these are now computed outside
if recordedPeak < 0
    recordedPeak = 0;
end

sig = sqrt(SIGMA(end,end));
D = muk(end)-recordedPeak;
alph = -D/sig;

% if D<=-6*sig
%     ecdi = 0;
% else
ecdi = D*(1-normcdf2(alph))+sig*normpdf(alph);
% end

% Stage Cost: energy charge (must be in kW)
alph2 = -muk(end)/sig;
ec = muk(end)*(1-normcdf2(alph2))+normpdf(alph2)*sig;

% Overall stage cost (includes damage)
stageCost = Rd*norm(sys.Dk(:,:,k)*uk) + Re*ec + peakCost*ecdi;

ns = size(xk,1);
% muk = [Ak*xk+Bk*uk+Fk*dk;mu];

%%%%%%%%%%%%% Done in another function now %%%%%%%%%%%%%%%%%%%%%%%%
% SIGMA = [   sys.Gk(:,:,k)*sys.Swk(:,:,k)*sys.Gk(:,:,k)'      sys.Gk(:,:,k)*sys.Swk(:,:,k)*sys.Hk(:,:,k)'
%             sys.Hk(:,:,k)*sys.Swk(:,:,k)*sys.Gk(:,:,k)'      sys.Hk(:,:,k)*sys.Swk(:,:,k)*sys.Hk(:,:,k)' + sys.Svk(:,:,k) ];
% SIGMA = SIGMA + tol*eye(3);
S11 = inv(SIGMA(1:ns,1:ns)-1/SIGMA(end,end)*SIGMA(1:ns,end)*SIGMA(1:ns,end)');
S12 = -1/SIGMA(end,end)*S11*SIGMA(1:ns,end);
S22 = 1/SIGMA(end,end)*(1+SIGMA(1:ns,end)'/SIGMA(end,end)*S11*SIGMA(1:ns,end));
SIGMAi = [S11,S12;S12',S22];
Hi1 = inv(GAMMAi+SIGMAi);
bias = GAMMAi*cm+SIGMAi*muk;
Gicm = GAMMAi*cm;

% dm = -sum(bias.*(Hi1*bias)) + sum(cm.*cm)/sig^2 + muk'*SIGMAi*muk;
vm = Hi1*bias;
dm = -sum(bias.*vm) + sum(cm.*Gicm) + muk'*SIGMAi*muk;
sigzz = sqrt(Hi1(end,end));
arg1 = (recordedPeak-vm(end,:))/sigzz;
I1 = exp(-0.5*dm).*(1-normcdf2(arg1));

% H = 1/sig^2*eye(ns)+S11;
H = GAMMAi(1:ns,1:ns)+S11;
Hi = inv(H);
b = S11*muk(1:ns);
c = 2*muk(1:ns)'*S12;
d = 2*S12'*Hi;
g3 = S22-S12'*Hi*S12;
% fm = (GAMMAi(1:ns,1:ns)*cm(1:ns,:)+b);
fm = Gicm(1:ns,:)+b; % \Gamma^{-1} c_{m,x} + S_{k,11}*\mu_{k,x}

g2 = d*fm-c;
g1 = GAMMAi(end)*(recordedPeak-cm(end,:)).^2 +...
    sum(cm(1:ns,:).*Gicm(1:ns,:)) + muk(1:ns)'*S11*muk(1:ns) -...
    sum(fm.*(Hi*fm));
arg2 = sqrt(g3)*(recordedPeak-muk(end)+g2/2/g3);
I2 = exp(-g1/2+g2.^2/8/g3).*normcdf2(arg2);

% bias = rbfVFA.theta(end,k+1);
EVkp1 = (sqrt(det(Hi1)/det(SIGMA))*I1+sqrt(1/g3/det(H)/det(SIGMA))*I2)*wm+wmb; % wmb=0 if bias is not present

ecost = stageCost + sys.gamma*EVkp1; % can't be negative