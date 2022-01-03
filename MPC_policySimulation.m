function [X,U,cost,Y,costComponents] = MPC_policySimulation(sys,sysb,dk,...
    trainModel,Nsamples,NH,t0,t1,W,V,plotFlag)

% unpack system
dt = sys.dt;
t = t0:dt:t1+dt;
N = length(t);
x0 = sys.x0;
quietMode = false;
na = size(trainModel.ar,1)-1;
Ns = sysb.Ns;
Ni = sysb.Ni;

% Generate noise sequences for sample average
% approximation. This is done once for the entire horizon
Wsa = nan(N,Nsamples);
Vsa = zeros(N,Nsamples);
for jj = 1:Nsamples
    DEMAND = simulateDemandAR(trainModel,t,x0(1));
    for kk = 1:N-1
        Ck = sys.Ck(:,:,kk);
        Ek = sys.Ek(:,:,kk);
        Hk = sys.Hk(:,:,kk);
        ybar = Ck(1)*DEMAND(kk)+Ek*dk(kk); % no control
        Wsa(kk,jj) = (DEMAND(kk+1) - ybar)/Hk(1);
    end
end

% Simulate
U = zeros(Ni,N-1);
cost = zeros(1,N-1); % stage costs
costComponents = zeros(3,1); % cost components
X = zeros(Ns+1,N);
Y = zeros(1,N-1);
X(:,1) = x0;
tic
dispInds = floor(linspace(1,N-1,10));
for ii=1:N-1
    % get state
    xk = X(1:end-1,ii);
    recordedPeak = X(end,ii);
    
    % finite horizon OCP 
    sysb.x0 = xk; % set initial condition to current state
    LEN = min([N-ii-1,NH]);
    D = dk(ii:ii+LEN);   
    [Yp,Xp,Up] = finiteHorizonOCP(sysb,ii,D,Wsa(ii:ii+LEN,:),...
        Vsa(ii:ii+LEN,:),recordedPeak);
    if ~quietMode && ismember(ii,dispInds)
        disp("k="+num2str(ii)+", Time: "+num2str(toc,'%.2f'))
    end
    uk = Up(1);
    U(:,ii) = uk;
    
    % apply input to system
    Ak = sys.Ak(:,:,ii);
    Bk = sys.Bk(:,:,ii);
    Ck = sys.Ck(:,:,ii);
    Dk = sys.Dk(:,:,ii);
    Ek = sys.Ek(:,:,ii);
    Fk = sys.Fk(:,:,ii);
    Gk = sys.Gk(:,:,ii);
    Hk = sys.Hk(:,:,ii);
    
    xpost = Ak*xk+Fk*dk(ii)+Bk*uk;
    yp = Ck*xk+Dk*uk+Ek*dk(ii);
    
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

