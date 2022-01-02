function [V,uk] = bellmanEquationLookup(k,xk,recordedPeak,dk,sys,lookup,cons)



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
uub = cons.uub;
ulb = cons.ulb;

% warning for case when disturbance affects bounds
ind = find(~all(Gk==0,2) & ~isinf(xlb) & ~isinf(-xlb) );
for jj = 1:length(ind)
    warning(['Bounds applied to state ',num2str(ind(jj)),...
        'will only be approximate since the disturbance',...
        'affects this state'])
end

% Set up constraints
% AA = [Bk;-Bk];
% bb = [xub-Ak*xk;Ak*xk-xlb];
% AA(isinf(bb),:) = [];
% bb(isinf(bb),:) = [];

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
obj = @(u) expectedCostLookup(k,u,muk,SIGMAk,recordedPeak,lookup,sys);
lbm = max([ulb;bb(AA<0)./AA(AA<0)]); % most restrictive LB
ubm = min([uub;bb(AA>0)./AA(AA>0)]); % most restrictive UB
% [~,V] = fminbnd(obj,lbm,ubm,opt2);
Nu = length(lookup.ugrid);
val = inf(1,Nu);
ug = lookup.ugrid;
for jj = 1:Nu
    if ug(jj)<=ubm && ug(jj)>=lbm
        val(jj) = obj(ug(jj));
    end
end
[V,idx] = min(val);
uk = lookup.ugrid(idx);
% V = max([0,V]);
%     [ecost,EVkp1,stageCost] = obj(uk);