function [Y,X,U] = finiteHorizonOCP(sysb,k0,D,W,V,recordedPeak)

K = size(W,1)-1; % horizon
Ns = sysb.Ns; % state dimension
Ni = sysb.Ni; % input dimension
Nd = sysb.Nd; % deterministic disturbance dimension
Nw = sysb.Nw; % noise dimension
Ntraj = size(W,2); % number of trajectories
x0 = reshape(sysb.x0,numel(sysb.x0),1);
x0 = x0(1:Ns); % sysb.x0 = x(k0) set outside

% extract system matrices from overall matrices
k1 = k0+K;
indOutput = k0:k1; % scalar output

% extract subsections of block matrices
% and append matrices to beginning to obtain x0
% as the first state

% dynamic equation matrices
indStates = (k0-1)*Ns+1:k1*Ns;
indInputs = (k0-1)*Ni+1:k1*Ni;
indDist = (k0-1)*Nd+1:k1*Nd;
indNoise = (k0-1)*Nw+1:k1*Nw;
Ab = sysb.Ab(indStates,:);
Ab = [eye(Ns);Ab(1:end-Ns,:)];
Bb = sysb.Bb(indStates,indInputs);
Fb = sysb.Fb(indStates,indDist);
Gb = sysb.Gb(indStates,indNoise);

% output matrices
Cb = sysb.Cb(indOutput,indStates);
Db = sysb.Db(indOutput,indInputs);
Eb = sysb.Eb(indOutput,indDist);
Hb = sysb.Hb(indOutput,indNoise);

BIG = 1e9;
% D = D(k0:k1);
b = Cb*Ab*x0 + (Cb*Fb+Eb)*D + (Cb*Gb+Hb)*W+ V;
b1 = Ab*x0+Fb*D+Gb*W;
XU = repmat(sysb.stateUpperLimits,[K+1,Ntraj]);
XL = repmat(sysb.stateLowerLimits,[K+1,Ntraj]);
UU = repmat(sysb.controlUpperLimits,[K+1,1]);
UL = repmat(sysb.controlLowerLimits,[K+1,1]);
Rd = sysb.damageCost*ones(1,K+1);
Re = sysb.energyCost*ones(1,K+1);

% process inf state constraints
XU = replaceInf(XU,BIG);
XL = replaceInf(XL,BIG);

% cvx_solver mosek
MAT = Cb*Bb+Db;
recordedPeak = repmat(recordedPeak,1,Ntraj);
% if sameSOC % has constraint that SOC(end)==SOC(start)
%     cvx_begin quiet
%     variable u(N+1)
%     minimize(sysb.peakCost/Ntraj*max([MAT*u+b;recordedPeak],[],2) +...
%         Rd*abs(Db*u) + 1/Ntraj*sum(Re*pos(MAT*u+b)))
%     subject to
%         UL <= u <= UU;
%         XL <= Bb*u+b1 <= XU;
%         Bb(end,:)*u + b1(end) >= sysb.x0(2);
%     cvx_end
% else
cvx_begin quiet
variable u(K+1)
minimize(sysb.peakCost/Ntraj*sum(max([repmat(MAT*u,1,Ntraj)+b;recordedPeak])) +...
    Rd*abs(Db*u) + ... % deterministic
    1/Ntraj*sum(sum(Re*pos(repmat(MAT*u,1,Ntraj)+b),1)))
subject to
    UL <= u <= UU;
    XL <= repmat(Bb*u,1,Ntraj)+b1 <= XU;
cvx_end
% end

U = u;
X = Ab*x0+Bb*U+Fb*D+Gb*W;
% XP = Ab*x0+Bb*U+Fb*D; % post decision state
% Y = Cb*X+Db*U+Eb*D+Hb*W+V;
Y = MAT*U+b;
Y = reshape(Y,[1,K+1,Ntraj]);
X = reshape(X,[Ns,K+1,Ntraj]);
% cost = cvx_optval;
