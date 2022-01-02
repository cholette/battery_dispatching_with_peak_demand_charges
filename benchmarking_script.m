% Comparing performance of RBF VFA, Lookup table, and Stochastic MPC on a
% full month simulation.
close all
delete(gcp('nocreate'))
clearvars
addpath demand_model utilities rbfVFA deterministic lookup_table
rng('shuffle')

% Number of workers for parfor loop. 
NW = 6; 
hpc = false; % true if you want to use the special HPC configuration below
saveResults = false; % save results as a .mat file? Overwritten to "true" if hpc == true

% Load in optimization parameters, RBFs and Lookup Tables. Please note the
% below are examples and may need to be changed for your configuration.
d = "./"; % base directory of saved files (use ./ for current directory)

% optimization model (obtained by running load_optimization_parameters.m)
opt_model = d+"optimization_model_january"; 

% array of strings for RBF VFA .mat files (each obtained by running backwardADP_top.m)
RBFFiles = d+["rbfBackwardADP_Jan_125_1"]; 

% array of strings for lookup table .mat files to simulate (each obtained by running 
% backwardDP_Lookup_top.m)
lookupFiles = d+["LookupBackwardDP_21_8000_1"];

% MPC parameters
Nsamples = 10; % Base number samples for MPC SAA, 10x this number is also simulated.
H = 48*2; % Number of 30-minute intervals for MPC optimisation horizon

% Setup CVX with Mosek. CVX for MATLAB is required for this script, but if
% you don't have Mosek, comment out the solver selection line. In our
% experience, Mosek is faster for this application than the default CVX
% solvers. The below code should be edited for your HPC setup, if you are
% using one.
if hpc % setup CVX with MOSEK
    cdir = pwd;
    cd ~/cvx/
    cvx_setup
    cvx_solver mosek
    cd(cdir)
    saveResults = true;
    plotResults = false;
	sname = "benchmarking_results.mat"
else
    cvx_clear
    cvx_solver mosek
    plotResults = true;
end

testMode = true; 
plotDuringSimulation = false;
verbosity = 'q'; % 'q' for quiet, anything else for verbose
N = 250; % Number of simulations to run
test = "random"; % "deterministic" or "random" initial state (initial peak is still zero)
x0 = nan(3,1); % x0 = [35;36;0]; % value only matters if test="deterministic".
t0 = datetime(2019,1,1);
file = 'data/Parameters_NO_AC.xlsx';

% load in demand data
disp("Loading "+opt_model)
dat = load(opt_model);
xub = dat.xub;
xlb = dat.xlb;
ulb = dat.ulb;
uub = dat.uub;
sys = dat.sys;
dk = dat.dd;
T = dat.T;
DAYS = dat.DAYS;
trainModel = dat.trainModel;
tracesTrain = dat.tracesTrain;
if testMode
    disp("Test Model. Using test trace from file: "+opt_model)
    testModel = dat.testModel;
    tracesTest = dat.tracesTest;
else
    disp("Evaluating on training model in file: "+opt_model)
end
sys.x0 = [];
sys.epsilon = 0; % greedy
t = dat.t;

% load in RBFs
Nrbf = length(RBFFiles);
rbfs = cell(1,Nrbf);
ylbs = []; yubs = []; lbs = []; ubs = [];
for ii = 1:Nrbf
    previousResultFile = RBFFiles(ii);
    disp("Initializing from " + previousResultFile)
    temp = load(previousResultFile);
    FF = fields(temp.rbfVFA);
    for jj = 1:length(FF)
        rbfVFA.(FF{jj}) = temp.rbfVFA.(FF{jj});
    end
    rbfs{ii} = rbfVFA;
    ylbs = [ylbs,temp.sys.outputLowerLimits];
    yubs = [yubs,temp.sys.outputUpperLimits];
    lb = [lbs,temp.lb];
    ub = [ubs,temp.ub];
    clear temp rbfVFA
end

t1 = t0 + hours(0.5*48*DAYS-0.5);

% load in lookup tables
Nlookup = length(lookupFiles);
lookupTables = cell(1,Nlookup);
for ii = 1:length(lookupFiles)
    dat = load(lookupFiles(ii));
    lookupTables{ii} = dat.lookupTable;
    clear dat
end

% initialize variables
costADPs = zeros(N,Nrbf); costPK = zeros(N,1); costUC = zeros(N,1); 
costMPC = zeros(N,1); costMPC2 = zeros(N,1); costLookup = zeros(N,Nlookup);

timeADPs = zeros(N,Nrbf); timePK = zeros(N,1); timeMPC = zeros(N,1);
timeMPC2 = zeros(N,1); timeLookup = zeros(N,Nlookup);

costComponentsADPs = zeros(N,Nrbf,3); costComponentsPK = zeros(N,1,3);
costComponentsUC = zeros(N,1,3); costComponentsMPC = zeros(N,1,3);
costComponentsMPC2 = zeros(N,1,3); costComponentsLookup = zeros(N,Nlookup,3);

peakADPs = zeros(N,Nrbf); peakPK = zeros(N,1); peakUC = zeros(N,1); 
peakMPC = zeros(N,1); peakMPC2 = zeros(N,1);peakLookup = zeros(N,Nlookup);

% make the stacked systems
sysb = makeStackedSystem(sys,T-1);

formatSpec = "n = %d";
formatSpec = formatSpec+"\n  x0: [%.1f, %.1f, %.1f]";
for ii = 1:Nrbf
    formatSpec = formatSpec+"\n  Controlled (ADP%d): [%.2f, %.1f]";
end
for ii = 1:Nlookup
    formatSpec = formatSpec+"\n  Controlled (Lookup%d): [%.2f, %.1f]";
end
formatSpec = formatSpec+"\n  MPC:               [%.2f, %.1f]";
formatSpec = formatSpec+"\n  MPC (10x):         [%.2f, %.1f]";
formatSpec = formatSpec+"\n  Perfect:           [%.2f, %.1f]";
formatSpec = formatSpec+"\n  Uncontrolled:      [%.2f, %.1f]\n";

% Handle parallel pool
if NW>1
    parpool(NW);
end

parfor nn = 1:N
    sysnn = sys;
    sysbnn = sysb;
    tstartLoop = tic;
    
    disp("Sampling initial state ... ")
    if test == "random"
        sysnn.x0 = sampleInitialState(sysnn,testModel,[xlb;0],[xub;min(ub(end,:))],true);
    elseif test == "deterministic"
        sysnn.x0 = x0;
    else
        error('test not recognized')
    end
    x0 = sysnn.x0; % same initial state for all simulations
    
    if testMode % generate demand profile from test model
        DEMAND = simulateDemandAR(testModel,[t,t(end)+duration(0,30,0)]);
        [~,E] = simulateDemandAR(testModel,t); % get expected value
    else % simulate on training model
        DEMAND = simulateDemandAR(trainModel,[t,t(end)+duration(0,30,0)]);
        [~,E] = simulateDemandAR(trainModel,t); % get expected value
    end
    
    DEMAND = DEMAND(:,1);
    W = nan(1,T);
    V = zeros(1,T);
    XADP = zeros(3,T+1,Nrbf);
    YADP = zeros(1,T,Nrbf);
    CADP = zeros(1,T+1,Nrbf);
    CCOMPS = zeros(3,Nrbf);
    XLK = zeros(3,T+1,Nlookup);
    YLK = zeros(1,T,Nlookup);
    CLK = zeros(1,T+1,Nlookup);
    CCOMPSLK = zeros(3,Nlookup);
    
    % Compute W and V sequences from demand (need to use training
    % model for a fair comparison, which must be PAR(1) for the below to work)
    for jj = 1:T
        Ck = sysnn.Ck(:,:,jj);
        Ek = sysnn.Ek(:,:,jj);
        Hk = sysnn.Hk(:,:,jj);
        ybar = Ck(1)*DEMAND(jj)+Ek*dk(jj); % no control
        W(jj) = (DEMAND(jj+1) - ybar)/Hk(1);
    end
    
    % Simualte RBF VFAs
    indPeaksADPs = nan(1,Nrbf);
    disp("Simulating ADP Policies ...")
    rbfsnn = rbfs;
    totTimeADP = 0;
    for ii=1:Nrbf
        tstart = tic;
        sysnnii = sysnn;
        disp("    Simulating RBF policy " + num2str(ii)+" (out of "+num2str(Nrbf)+")")
        sysnnii.outputUpperLimits = yubs(:,ii);
        sysnnii.outputLowerLimits = ylbs(:,ii);       
        [XADP(:,:,ii),~,CADP(1,:,ii),YADP(1,:,ii),~,~,CCOMPS(:,ii)] = ...
            simulatePolicy_RBFVFA(sysnnii,dk,rbfsnn{ii},t0,t1,[],verbosity,W,V,plotDuringSimulation);
        timeADPs(nn,ii) = toc(tstart);
        totTimeADP = totTimeADP + timeADPs(nn,ii);
    end
    disp("... done. Average ADP time: "+num2str(totTimeADP/Nrbf/60)+" mins ")
    
    % Compute costs and peaks for ADP
    for ii = 1:Nrbf
        costComponentsADPs(nn,ii,:) = CCOMPS(:,ii);
        costADPs(nn,ii) = sum(CADP(1,:,ii),2);
        [p,ind] = max(YADP(1,:,ii));
        peakADPs(nn,ii) = max([sysnn.x0(end),p]);
        indPeaksADPs(ii) = ind;
    end
    
    % Compute Perfect control
    disp("Calculating Perfect Knowledge Control...")
    D = reshape(dk,numel(dk),1);
    sysbnn.x0 = sysnn.x0(1:2);
    tstart = tic;
    [XPK,~,UPK,COSTPK,YPK,stageCostPK,CCOMPSPK] = simulatePerfectControl(...
        sysbnn,D,W',V',sysnn.x0(end),false);
    timePK(nn) = toc(tstart);
    costComponentsPK(nn,1,:) = CCOMPSPK;
    disp("... done. "+num2str(timePK(nn)/60)+" mins ")
    
    % Simulate Lookups
    indPeaksLookup = nan(1,Nlookup);
    disp("Simulating Lookup 1 ...")
    lookupTablesnn = lookupTables;
    totTimeLookup = 0;
    for mm=1:Nlookup
        tstart = tic;
        sysnnii = sysnn;
        disp("    Simulating lookup policy " + num2str(mm)+" (out of "+num2str(Nlookup)+")")
        sysnnii.outputUpperLimits = yubs(:,mm);
        sysnnii.outputLowerLimits = ylbs(:,mm);       
        [XLK(:,:,mm),~,CLK(1,:,mm),YLK(1,:,mm),~,~,CCOMPSLK(:,mm)] = ...
            simulatePolicyLookupPAR1(sysnnii,dk,lookupTablesnn{mm},t0,t1,...
            verbosity,W,V,plotDuringSimulation);
        temp = toc(tstart);
        timeLookup(nn,mm) = temp;
        totTimeLookup = totTimeLookup + temp;
    end
    disp("... done. Average Lookup Table time: "+num2str(totTimeLookup/Nlookup/60)+" mins ")
    
    % Costs and peaks for Lookup
    for mm = 1:Nlookup
        costComponentsLookup(nn,mm,:) = CCOMPSLK(:,mm);
        costLookup(nn,mm) = sum(CLK(1,:,mm),2);
        [p,ind] = max(YLK(1,:,mm));
        peakLookup(nn,mm) = max([x0(end),p]);
        indPeaksLookup(mm) = ind;
    end
    
      
    % Simualte Stochastic MPC
    disp("Simulating MPC ...")
    tstart = tic;
    [XMPC,~,COSTMPC,YMPC,CCOMPSMPC] = MPC_policySimulation(sysnn,sysbnn,dk',...
        trainModel,Nsamples,H,t0,t1,W,V,plotDuringSimulation);
    timeMPC(nn) = toc(tstart);
    costComponentsMPC(nn,1,:) = CCOMPSMPC;
    disp("... done. "+num2str(timeMPC(nn)/60)+" mins ")
    
    % Simulate Stochast MPC with 10x samples for SAA
    disp("Simulating MPC(10x) ...")
    tstart = tic;
    [XMPC2,~,COSTMPC2,YMPC2,CCOMPSMPC2] = MPC_policySimulation(sysnn,sysbnn,dk',...
        trainModel,10*Nsamples,H,t0,t1,W,V,plotDuringSimulation);
    timeMPC2(nn) = toc(tstart);
    costComponentsMPC2(nn,1,:) = CCOMPSMPC2;
    disp("... done. "+num2str(timeMPC2(nn)/60)+" mins ")
     
    % costs of PK and MPC
    costPK(nn) = COSTPK;
    costMPC(nn) = sum(COSTMPC);
    costMPC2(nn) = sum(COSTMPC2);
    
    % uncontrolled, perfect, and MPC
    [ucdMax1,indUC] = max(XADP(1,2:end,1)); % should be the same for all simulations (check)
    [cdMaxPK1,indPK] = max(YPK);
    [cdMaxMPC1,indMPC] = max(YMPC);
    [cdMaxMPC12,indMPC2] = max(YMPC2);
    peakUC(nn) = max([sysnn.x0(end),ucdMax1]);
    peakPK(nn) = max([sysnn.x0(end),cdMaxPK1]);
    peakMPC(nn) = max([sysnn.x0(end),cdMaxMPC1]);
    peakMPC2(nn) = max([sysnn.x0(end),cdMaxMPC12]);
    
    % uncontrolled
    YUC = XADP(1,2:end,1);
    costComponentsUC(nn,:) = [dot(sys.energyCost*(YUC>0),YUC),0,...
        peakUC(nn)*sysnn.peakCost]; % [energy cost,damage cost, peak cost]
    costUC(nn) = dot(sys.energyCost*(YUC>0),YUC) + ...
        peakUC(nn)*sysnn.peakCost; % only works if gamma is one!
    
    %%%%%%%%%%%%%% Display %%%%%%%%%%%%%%%%%%%
    A = nn;
    A = [A,reshape(x0,1,3)];
    for ii = 1:Nrbf
        A = [A,ii,peakADPs(nn,ii),costADPs(nn,ii)];
    end

    for mm = 1:Nlookup
        A = [A,mm,peakLookup(nn,mm),costLookup(nn,mm)];
    end
    A = [A,peakMPC(nn),costMPC(nn),peakMPC2(nn),costMPC2(nn),...
        peakPK(nn),costPK(nn),peakUC(nn),costUC(nn)];
    fprintf(formatSpec,A)
    disp("  Loop time: "+num2str(toc(tstartLoop))+" secs")
    
    if plotResults
        % SOC
        c = lines(Nrbf+Nlookup+4);
        figure(1)
        for ii = 1:Nrbf
            plot(t,XADP(2,2:end,ii),'Color',c(ii+1,:),"DisplayName","ADP"+num2str(ii))
            hold on
        end
        for mm = 1:Nlookup
            plot(t,XLK(2,2:end,mm),'Color',c(mm+Nrbf,:),"DisplayName","Lookup"+num2str(mm))
            hold on
        end
        plot(t,XPK(2,1:end),':','Color',c(end-2,:),"DisplayName","Perfect");
        plot(t,XMPC(2,2:end),'--','Color',c(end-1,:),"DisplayName","MPC"); 
        plot(t,XMPC2(2,2:end),'--','Color',c(end,:),"DisplayName","MPC(10x)"); hold off
        ylabel('Battery stored energy (kWh)')
        legend('location','Southeast')
        hold off
        drawnow
        
        % Demand
        figure(2)
        plot(t,YUC,'Color',c(1,:),"DisplayName","Uncontrolled"); hold on;
        plot(t(indUC),ucdMax1,'*','Color',c(1,:),"DisplayName","UC Peak")
        Dlast = XADP(1,:,1);
        for ii = 1:Nrbf
            % check that demand is exactly the same
            if ~all(XADP(1,:,ii) == Dlast)
                error('Demand not the same')
            end
            Dlast = XADP(1,:,ii);
            
            plot(t,YADP(1,:,ii),'Color',c(ii+1,:),"DisplayName","ADP"+num2str(ii))
            plot(t(indPeaksADPs(ii)),peakADPs(nn,ii),"*",'Color',...
                c(ii+1,:),"DisplayName","ADP"+num2str(ii)+" Peak")
        end
        for mm = 1:Nlookup
            % check that demand is exactly the same
            if ~all(XLK(1,:,mm) == Dlast)
                error('Demand not the same')
            end
            Dlast = XLK(1,:,mm);
            
            plot(t,YLK(1,:,mm),'Color',c(mm+Nrbf,:),"DisplayName","Lookup"+num2str(mm))
            plot(t(indPeaksLookup(mm)),peakLookup(nn,mm),"*",'Color',...
                c(mm+1,:),"DisplayName","Lookup"+num2str(mm)+" Peak")
        end
        plot(t,YPK,':','Color',c(end-2,:),"DisplayName","PK")
        plot(t(indPK),cdMaxPK1,'*','Color',c(end-2,:),"DisplayName","PK Peak")
        plot(t,YMPC,'--','Color',c(end-1,:),"DisplayName","MPC")
        plot(t(indMPC),cdMaxMPC1,'*','Color',c(end-1,:),"DisplayName","MPC Peak")
        plot(t,YMPC2,'--','Color',c(end,:),"DisplayName","MPC(10x)")
        plot(t(indMPC2),cdMaxMPC12,'*','Color',c(end,:),"DisplayName","MPC Peak(10x)")
        legend('location','Southeast')
        ylabel('Demand (kW)'); hold off
        hold off
        drawnow
        
        % Cumulative maximums
        figure(3)
        plot(t,cummax(YUC),'Color',c(1,:),"DisplayName","Uncontrolled"); hold on;
        for ii = 1:Nrbf
            plot(t,cummax(YADP(1,:,ii)),'Color',c(ii+1,:),"DisplayName","ADP"+num2str(ii))
        end
        for mm = 1:Nlookup
            plot(t,cummax(YLK(1,:,mm)),'Color',c(mm+Nrbf,:),"DisplayName","Lookup"+num2str(mm))
        end
        plot(t,cummax(YPK),':','Color',c(end-2,:),"DisplayName","PK")
        plot(t,cummax(YMPC),'--','Color',c(end-1,:),"DisplayName","MPC")
        plot(t,cummax(YMPC2),'--','Color',c(end,:),"DisplayName","MPC(10x)")
        legend('location','Southeast')
        ylabel('Cumulative Peak Demand (kW)'); hold off
        hold off
        drawnow
    end
end

%  summary statistics
formatSpec = split(formatSpec,"\n");
formatSpec(1:2) = [];
formatSpec = join(formatSpec,"\n");
A=[];
for ii = 1:Nrbf
    A = [A,ii,mean(peakADPs(:,ii)),mean(costADPs(:,ii))];
end
for mm = 1:Nlookup
    A = [A,mm,mean(peakLookup(:,mm)),mean(costLookup(:,mm))];
end
A = [A,mean(peakMPC),mean(costMPC),mean(peakMPC2),mean(costMPC2),...
    mean(peakPK),mean(costPK),mean(peakUC),mean(costUC)];
disp('=================== Summary (Means) =========================')
fprintf(formatSpec,A)

if saveResults
    save(sname,'-v7.3')
    disp("File saved as: "+sname)
end

delete(gcp('nocreate'))