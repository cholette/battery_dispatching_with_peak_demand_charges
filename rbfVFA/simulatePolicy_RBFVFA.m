function [X,U,cost,Y,W,V,costComponents] = simulatePolicy_RBFVFA(sys,dk,rbfVFA,t0,t1,varargin)

% unpack system
dt = sys.dt;
t = t0:dt:t1;
N = length(t);
x0 = sys.x0;
quietMode = false;
% ug = linspace(-1,1,100);

if nargin > 5 && ~isempty(varargin{1})
    Vf = varargin{1};
else
    Vf = @(x) 0*x(end);
end

if nargin > 6
    varargin{2} = char(varargin{2}(1));
    if varargin{2}(1)=="q" % not in quiet mode
        quietMode = true;
    end
end

% if trace is supplied, use it as the true demand (instead of demand model
% embedded in sys)
if nargin > 7 % test trace supplied
    if ~quietMode
        disp("Testing ADP on a specified demand trace...")
    end
    testMode = true;
    W = varargin{3};
    V = varargin{4};
else
    if ~quietMode 
        disp("Using same system in optimization and simulation")
    end
    testMode = false;
    W = nan(N,1);
    V = nan(N,1);
end

if nargin > 9
    plotFlag = varargin{5};
    if plotFlag 
        h = figure("Name","Lookup Policy");
    end
else
    plotFlag = false;
end

Ns = size(sys.Ak(:,:,1),2);
Ni = size(sys.Bk(:,:,1),2);
% Nu = length(ug);

% constraints
cons.xub = sys.stateUpperLimits;
cons.xlb = sys.stateLowerLimits;
cons.uub = sys.controlUpperLimits;
cons.ulb = sys.controlLowerLimits;
cons.yub = sys.outputUpperLimits;
cons.ylb = sys.outputLowerLimits ;

% Simulate
U = zeros(Ni,N);
cost = zeros(1,N+1); % stage costs
costComponents = zeros(3,1); % cost components
X = zeros(Ns+1,N+1);
Y = zeros(1,N);
X(:,1) = x0;
opt2 = optimset('display','off','tolX',1e-5,'FunValCheck','off');
for ii=1:N
    
    % get system matrices for optimization
    Ak = sys.Ak(:,:,ii);
    Bk = sys.Bk(:,:,ii);
    Ck = sys.Ck(:,:,ii);
    Dk = sys.Dk(:,:,ii);
    Ek = sys.Ek(:,:,ii);
    Fk = sys.Fk(:,:,ii);
    Gk = sys.Gk(:,:,ii);
    Hk = sys.Hk(:,:,ii);
       
    % get state
    xk = X(1:end-1,ii);
    recordedPeak = X(end,ii);
    
    if sys.epsilon > 0 && rand <= sys.epsilon % random
        
        % State constraints
        AAx = [Bk;-Bk];
        bbx = [cons.xub-Ak*xk-Fk*dk(ii);Ak*xk+Fk*dk(ii)-cons.xlb];
        AAx(isinf(bbx),:) = [];
        bbx(isinf(bbx),:) = [];

        % Output constraints
        AAy = [Dk;-Dk];
        bby = [cons.yub-Ck*xk-Ek*dk(ii);Ck*xk+Ek*dk(ii)-cons.ylb];
        AAy(isinf(bby),:) = [];
        bby(isinf(bby),:) = [];
        AA = [AAx;AAy];
        bb = [bbx;bby];
        
        uk = cons.ulb + (cons.uub-cons.ulb)*rand;
%         uk = (uub+ulb)/2 + (uub-ulb)/6*randn;
        uk = 1*(uk>1) + (uk<=1)*uk;
        uk = -1*(uk<-1) + (uk>=-1)*uk;
        if any( AA*uk >= bb ) % find closest point that satisfies the contraints
            obj = @(u) sum((u-uk).^2);
            lbm = max([cons.ulb;bb(AA<0)./AA(AA<0)]); % most restrictive LB
            ubm = min([cons.uub;bb(AA>0)./AA(AA>0)]); % most restrictive UB
            uk = fminbnd(obj,lbm,ubm,opt2);
            %             uk = fmincon(obj,0,AA,bb,[],[],ulb,uub,[],opt);
        end
    else % greedy
       [~,uk] = bellmanEquation_RBFVFA(ii,xk,recordedPeak,dk(ii),sys,rbfVFA,cons);
    end
    U(:,ii) = uk;
    
    % apply input to system
    xpost = Ak*xk+Fk*dk(ii)+Bk*uk;
    yp = Ck*xk+Dk*uk+Ek*dk(ii);
    
    % get system matrices for simulation
    if ~testMode  % generate noise from model statistics
        W(ii) = sqrt(sys.Swk(:,:,ii))*randn;
        V(ii) = sqrt(sys.Svk(:,:,ii))*randn;
    end
    
    % get process and measurement noise
    xp = xpost+Gk*W(ii);
    
    % get demand w control
    Y(ii) = yp+Hk*W(ii)+ V(ii);
    
    X(:,ii+1) = [   xp
        max([X(end,ii),Y(ii)])];
    
    % stage cost
    costComponents(1) = costComponents(1) + sys.gamma^(ii-1)*sys.energyCost*(Y(ii)>0)*Y(ii); % energy cost
    costComponents(2) = costComponents(2) + sys.gamma^(ii-1)*sys.damageCost*norm(sys.Dk(:,:,ii)*uk);
    costComponents(3) = costComponents(3) + sys.gamma^(ii-1)*sys.peakCost*( (Y(ii)>X(end,ii))*(Y(ii)- recordedPeak) );
    cost(ii) = sys.gamma^(ii-1)*( sys.damageCost*norm(sys.Dk(:,:,ii)*uk) +...
        sys.energyCost*(Y(ii)>0)*Y(ii) + ...
        sys.peakCost*( (Y(ii)>X(end,ii))*(Y(ii)- recordedPeak) ));
    
    if plotFlag
       tiledlayout(4,1)
       
       nexttile       
       plot(t(1:ii-1),X(1,2:ii),"DisplayName","Uncontrolled")
       hold on
       plot(t(1:ii),Y(1:ii),"DisplayName","Controlled")
       title("Net Demand")
       legend("location","NortheastOutside")
       
       nexttile
       plot(t(1:ii),U(1:ii))
       title("Control")
       
       nexttile
       plot(t(1:ii-1),X(2,2:ii))
       title("SOC")
       
       nexttile
       plot(t(1:ii-1),X(3,2:ii))
       title("Peak Demand")
       
       drawnow
    end
end

cost(N+1) = Vf(X(:,N+1));
