function [X,U,cost,Y,W,V,costComponents] = simulatePolicyLookupPAR1(sys,dk,lookup,t0,t1,varargin)

% unpack system
dt = sys.dt;
t = t0:dt:t1;
N = length(t);
x0 = sys.x0;
% d = length(lookupVFA.ind);

quietMode = false;
if nargin > 5
    varargin{1} = char(varargin{1}(1));
    if varargin{1}(1)=="q" % quiet mode
        quietMode = true;
    end
end

% if trace is supplied, use it as the true demand (instead of demand model
% embedded in sys)
if nargin > 6 % test trace supplied
    if ~quietMode
        disp("Testing ADP on a specified demand trace...")
    end
    testMode = true;
    W = varargin{2};
    V = varargin{3};
else
    if ~quietMode 
        disp("Using same system in optimization and simulation")
    end
    testMode = false;
    W = nan(N,1);
    V = nan(N,1);
end

if nargin > 7
    plotFlag = varargin{4};
    if plotFlag 
        h = figure("Name","Lookup Policy");
    end
else
    plotFlag = false;
end

Ns = size(sys.Ak(:,:,1),2);
Ni = size(sys.Bk(:,:,1),2);

% constraints
cons.xub = sys.stateUpperLimits;
cons.xlb = sys.stateLowerLimits;
cons.uub = sys.controlUpperLimits;
cons.ulb = sys.controlLowerLimits;
cons.yub = inf;
cons.ylb = -inf;
% cons.yub = sys.outputUpperLimits;
% cons.ylb = sys.outputLowerLimits ;

% Simulate
U = zeros(Ni,N);
cost = zeros(1,N+1); % stage costs
costComponents = zeros(3,1); % cost components
X = zeros(Ns+1,N+1);
Y = zeros(1,N);
% XP = zeros(d,N+1);
X(:,1) = x0;
for ii=1:N
    
    % get system matrices
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
    
    % get input by optimising value function
    [~,uk] = bellmanEquationLookupPAR1(ii,xk,recordedPeak,dk(ii),sys,lookup,cons);
    
    U(:,ii) = uk;
    
    % apply input to system
    xpost = Ak*xk+Fk*dk(ii)+Bk*uk;
    if ~testMode
        W(ii) = sqrt(sys.Swk(:,:,ii))*randn;
        V(ii) = sqrt(sys.Svk(:,:,ii))*randn;
    end
    xp = xpost+Gk*W(ii);
    
    % get demand w control
    yp = Ck*xk+Dk*uk+Ek*dk(ii);
    Y(ii) = yp+Hk*W(ii)+V(ii);
    
    % next state
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
       plot(h,t(1:ii-1),X(1,2:ii),"DisplayName","Uncontrolled")
       hold on
       plot(h,t(1:ii),Y(1:ii),"DisplayName","Controlled")
       title("Net Demand")
       legend("location","NortheastOutside")
       
       nexttile
       plot(h,t(1:ii),U(1:ii))
       title("Control")
       
       nexttile
       plot(h,t(1:ii-1),X(2,2:ii))
       title("SOC")
       
       nexttile
       plot(h,t(1:ii-1),X(3,2:ii))
       title("Peak Demand")
       
       drawnow
    end
    
    % stage cost
%     cost(ii) = sys.gamma^(ii-1)*( sys.damageCost*norm(sys.Dk(:,:,ii)*uk) +...
%         sys.energyCost*(Y(ii)>0)*Y(ii) + ...
%         sys.peakCost*( (Y(ii)>X(end,ii))*(Y(ii)- recordedPeak) ));
end