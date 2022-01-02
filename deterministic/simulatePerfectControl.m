function [X,XP,U,cost,Y,stageCost,costComponents] = simulatePerfectControl(sysb,...
    D,W,V,recordedPeak,sameSOC)

N = length(W)-1;
Ns = size(sysb.Ab,2);
x0 = reshape(sysb.x0,numel(sysb.x0),1);
x0 = x0(1:Ns);

Ab = sysb.Ab;
Bb = sysb.Bb;
Cb = sysb.Cb;
Db = sysb.Db;
Eb = sysb.Eb;
Fb = sysb.Fb;
Gb = sysb.Gb;
Hb = sysb.Hb;

BIG = 1e9;
b = Cb*Ab*x0 + (Cb*Fb+Eb)*D + (Cb*Gb+Hb)*W+ V;
b1 = Ab*x0+Fb*D+Gb*W;
XU = repmat(sysb.stateUpperLimits,[N+1,1]);
XL = repmat(sysb.stateLowerLimits,[N+1,1]);
UU = repmat(sysb.controlUpperLimits,[N+1,1]);
UL = repmat(sysb.controlLowerLimits,[N+1,1]);
% % Rb = makeBlock(sysb.alph,N+1);
Rd = sysb.damageCost*ones(1,N+1);
Re = sysb.energyCost*ones(1,N+1);

% process inf state constraints
XU = replaceInf(XU,BIG);
XL = replaceInf(XL,BIG);

% cvx_solver mosek
if sameSOC % has constraint that SOC(end)==SOC(start)
    cvx_begin quiet
    variable u(N+1)
    minimize(sysb.peakCost*norm([(Cb*Bb+Db)*u+b;recordedPeak],inf) +...
        Rd*abs(Db*u) + Re*pos((Cb*Bb+Db)*u+b))
    subject to
        UL <= u <= UU;
        XL <= Bb*u+b1 <= XU;
        Bb(end,:)*u + b1(end) >= sysb.x0(2);
    cvx_end
else
    cvx_begin quiet
    variable u(N+1)
    minimize(sysb.peakCost*norm([(Cb*Bb+Db)*u+b;recordedPeak],inf) +...
        Rd*abs(Db*u) + Re*pos((Cb*Bb+Db)*u+b))
    subject to
        UL <= u <= UU;
        XL <= Bb*u+b1 <= XU;
    cvx_end
end

U = u;
X = sysb.Ab*x0+sysb.Bb*U+Fb*D+Gb*W;
XP = Ab*x0+Bb*U+Fb*D; % post decision state
Y = Cb*X+Db*U+Eb*D+Hb*W+V;
X = reshape(X,[Ns,N+1]);
XP = reshape(XP,[Ns,N+1]);
cost = cvx_optval;

% compute stage costs
% cost(ii) = sys.gamma^(ii-1)*( sys.damageCost*norm(uk) + sys.energyCost*pos(uk) + ...
%         sys.peakCost* ( max( [Y(ii),X(end,ii)]) - X(end,ii) ));

M = cummax([recordedPeak;Y]);
XP = [XP;M(1:end-1)'];
costComponents(1) = Re*pos(Y); % energy cost
costComponents(2) = Rd*abs(Db*U); % damage
costComponents(3) = sysb.peakCost*norm([Y;recordedPeak],inf); % peak cost
stageCost = (sysb.peakCost*diff(M)+Rd'.*abs(Db*U) + Re'.*pos(Y))';
