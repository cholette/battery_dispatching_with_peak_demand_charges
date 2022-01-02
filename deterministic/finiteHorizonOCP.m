function [Y,X,U] = finiteHorizonOCP(sysb,k0,D,W,V,recordedPeak)

K = size(W,1)-1; % horizon
Ns = sysb.Ns; % state dimension
Ni = sysb.Ni; % input dimension
Nd = sysb.Nd; % deterministic disturbance dimension
Nw = sysb.Nw; % noise dimension
Ntraj = size(W,2); % number of trajectories
x0 = reshape(sysb.x0,numel(sysb.x0),1);
x0 = x0(1:Ns);

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
% Bb = [zeros(Ns,size(Bb,2));Bb(1:end-Ns,:)];
% Bb(:,end) = zeros(size(Bb,1),Ni);

Fb = sysb.Fb(indStates,indDist);
% Fb = [zeros(Ns,size(Fb,2));Fb(1:end-Ns,:)];
% Fb(:,end) = zeros(size(Fb,1),Nd);
Gb = sysb.Gb(indStates,indNoise);
% Gb = [zeros(Ns,size(Gb,2));Gb(1:end-Ns,:)];
% Gb(:,end) = zeros(size(Gb,1),Ni);

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
% % Rb = makeBlock(sysb.alph,N+1);
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

% XP = reshape(XP,[Ns,N+1,Ntraj]);
% cost = cvx_optval;

% compute stage costs
% cost(ii) = sys.gamma^(ii-1)*( sys.damageCost*norm(uk) + sys.energyCost*pos(uk) + ...
%         sys.peakCost* ( max( [Y(ii),X(end,ii)]) - X(end,ii) ));

% M = cummax([recordedPeak*ones(1,Ntraj);Y]);
% % XP = [XP;M(1:end-1)'];
% costComponents(1) = Re/Ntraj*sum(pos(Y),'all'); % average energy cost
% costComponents(2) = Rd*abs(Db*U); % damage
% costComponents(3) = sysb.peakCost/Ntraj*sum(max([repmat(MAT*u,1,Ntraj)+b;repmat(recordedPeak,1,Ntraj)],[],2)); % sysb.peakCost*norm([Y;recordedPeak],inf); % peak cost
% stageCost = (sysb.peakCost*diff(M)+Rd'.*abs(Db*U) + Re'.*pos(Y))';
