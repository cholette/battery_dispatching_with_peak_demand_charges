function [V,u] = bellmanEquation_RBFVFA(k,xk,recordedPeak,dk,sys,rbfVFA,cons,varargin)

opt2 = optimset('display','off','tolX',1e-7,'FunValCheck','off');
Nu = 101;

% Get system matrices & state
Ak = sys.Ak(:,:,k);
Bk = sys.Bk(:,:,k);
Ck = sys.Ck(:,:,k);
Gk = sys.Gk(:,:,k);
Fk = sys.Fk(:,:,k);
Ek = sys.Ek(:,:,k);
Dk = sys.Dk(:,:,k);
Hk = sys.Hk(:,:,k);

% unpack constraints
xub = cons.xub;
xlb = cons.xlb;
yub = cons.yub;
ylb = cons.ylb;
uub = cons.uub;
ulb = cons.ulb;

if nargin > 7 && ~isempty(varargin{1}) && varargin{1}(1)~='q'
    % warning for case when disturbance affects bounds
    ind = find(~all(Gk==0,2) & ~isinf(xlb) & ~isinf(-xlb) );
    for jj = 1:length(ind)
        warning(['Bounds applied to state ',num2str(ind(jj)),...
            'will only be approximate since the disturbance',...
            'affects this state'])
    end
    
    if ~all(Hk==0,"all") || ~all(sys.Svk(:,:,ii)==0,"all")
        warning(['Bounds applied to output will only be approximate '...
            'since there are (random) disturbances that affect it.'])
    end
    
end

% State constraints
AAx = [Bk;-Bk];
bbx = [cons.xub-Ak*xk-Fk*dk;Ak*xk+Fk*dk-cons.xlb];
AAx(isinf(bbx),:) = [];
bbx(isinf(bbx),:) = [];

% Output constraints
AAy = [Dk;-Dk];
bby = [cons.yub-Ck*xk-Ek*dk;Ck*xk+Ek*dk-cons.ylb];
AAy(isinf(bby),:) = [];
bby(isinf(bby),:) = [];
AA = [AAx;AAy];
bb = [bbx;bby];


% Solve Bellman Equation and save optimal value function
[muk,SIGMAk] = extractSystemInfo(sys,xk,dk,k);
obj = @(u) expectedCost_RBFVFA(k,u,xk,recordedPeak,sys,rbfVFA,muk,SIGMAk);
lbm = max([ulb;bb(AA<0)./AA(AA<0)]); % most restrictive LB
ubm = min([uub;bb(AA>0)./AA(AA>0)]); % most restrictive UB
if lbm >= ubm
    error("no feasible action")
    lb
    ub
end

[u,V] = fminbnd(obj,lbm,ubm,opt2);
V = max([0,V]);
