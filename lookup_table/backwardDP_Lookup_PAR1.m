function backwardDP_Lookup_PAR1(ncs,nw,plotFlag,NU,optModel,varargin)

if nargin>5
    restartFromFile = varargin{1};
    fileName = varargin{2};
else
    restartFromFile = false;
end

if ~restartFromFile
    rng('shuffle')
    delete(gcp('nocreate'))
    parpool(nw)
    
    %% Setup up user options
    ub = [130;216;130]; % Demand,Ebatt,RecordedPeak
    lb = [-60;36;0];
    nc = prod(ncs);
    sname = "LookupBackwardDP_"+num2str(NU)+"_"+num2str(nc)+".mat";
    
    %% Cost structure
    disp("LB:")
    lb
    disp("UB:")
    ub
    
	disp("Loading "+optModel)
	dat = load(optModel);
    xub = dat.xub;
    xlb = dat.xlb;
    ulb = dat.ulb;
    uub = dat.uub;
    ugrid = linspace(ulb,uub,NU); % u grid
    sys = dat.sys;
    dd = dat.dd;
    T = dat.T;
    DAYS = dat.DAYS;
    t0 = dat.t0;
    t1 = dat.t1;
    t = t0:hours(0.5):t1;
    trainModel = dat.trainModel;
    disp("Energy cost: "+num2str(sys.energyCost)+", Peak Cost: "+ ...
        num2str(sys.peakCost)+", Damage Cost: "+num2str(sys.damageCost));
    
    %% Set up Lookup table
    lookupTable = createLookupTable(ncs,lb,ub,T+1);
    Xs = lookupTable.xi(:,1:end-1)';
    Zs = lookupTable.xi(:,end);
    kStart = T;
else
    load(fileName)
    kStart = k-1;
    disp("Restarting from file "+fileName)
end

% Backward Induction
saveCounter = 0;
for k=kStart:-1:1 % time loop
    disp("k = "+num2str(k)+"...")
    tic
    
    % Get system matrices & state
    Ak = sys.Ak(:,:,k);
    Bk = sys.Bk(:,:,k);
    Ck = sys.Ck(:,:,k);
    Dk = sys.Dk(:,:,k);
    Ek = sys.Ek(:,:,k);
    Fk = sys.Fk(:,:,k);
    Gk = sys.Gk(:,:,k);
    Hk = sys.Hk(:,:,k);
    dk = dd(k);
    
    vv = nan(lookupTable.Ng,1);
    parfor ii = 1:lookupTable.Ng % Loop over samples
        
        xk = Xs(:,ii); %Xs(1:end-1,ii);
        recordedPeak = Zs(ii); %Xs(end,ii);
        
        % Set up constraints
        AA = [Bk;-Bk];
        bb = [xub-Ak*xk;Ak*xk-xlb];
        AA(isinf(bb),:) = [];
        bb(isinf(bb),:) = [];
        
        % Solve Bellman Equation and save optimal value
        [muk,S] = extractSystemInfo(sys,xk,dk,k);
        obj = @(u) expectedCostLookupPAR1(k,u,muk,S,recordedPeak,lookupTable,sys);
        lbm = max([ulb;bb(AA<0)./AA(AA<0)]); % most restrictive LB
        ubm = min([uub;bb(AA>0)./AA(AA>0)]); % most restrictive UB
        val = inf(1,NU);
        for jj = 1:NU
            if ugrid(jj)<=ubm && ugrid(jj)>=lbm
                val(jj) = obj(ugrid(jj));
            end
        end
        [v,idx] = min(val);
        %         uk = ugrid(idx);
        vv(ii) = v;
    end
    lookupTable.values(:,k) = vv;
    disp(" ... done! " + num2str(toc,'%.2f') + " seconds")
    
    if plotFlag
        scatter3(lookupTable.xi(:,1),lookupTable.xi(:,2),...
            lookupTable.xi(:,3),[],lookupTable.values(:,k))
        xlabel("Last demand")
        ylabel("SOC")
        zlabel("Recorded Peak")
        c = colorbar;
        c.Label.String = 'Cost-to-go';
        drawnow
    end
    
    saveCounter = saveCounter +1;
    if saveCounter/T>0.05
        save("backwardInductionState.mat_"+num2str(prod(ncs)))
        saveCounter = 0;
    end
end

%% TEST
% DAYS = 1;
t0 = datetime(2019,1,1);
t1 = t0 + hours(T*0.5-0.5);
t = t0:hours(0.5):t1;
N = 100;
x0=[35;38;0];

% datADP.sys = sys;
% datADP.sys.epsilon = 0; % greedy
% t1 = t0 + hours(0.5*48*DAYS-0.5);

% simulate
COST = zeros(1,N);
[~,indg] = aggragationFunction(x0',lookupTable);
PRED = lookupTable.values(indg);
sys.x0 = x0;
lookupTable.ugrid = ugrid;
disp('simulating...')
parfor nn = 1:N
    nn
    [X1,U1,COST1,Y1] = simulatePolicyLookupPAR1(sys,dd,lookupTable,t0,t1);
    COST(nn) = sum(COST1);
end

disp('========= Summary ===========================')
UB = mean(COST)+1.96/sqrt(N)*std(COST);
LB = mean(COST)-1.96/sqrt(N)*std(COST);
disp(['Controlled Cost: [',num2str(LB,'%.2f'),', ',num2str(UB,'%.2f'),']'])
disp(['Predicted Cost (DP): ',num2str(PRED,'%.2f')])

save(sname)
disp("Model saved as "+sname)
delete(gcp) % shut down parallel pool
