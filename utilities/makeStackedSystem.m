function sysb = makeStackedSystem(sys,N)

% Output equation first
Cb = makeBlock(sys.Ck);
Db = makeBlock(sys.Dk);
Eb = makeBlock(sys.Ek);
Hb = makeBlock(sys.Hk);

% now for the state equation
Ns = size(sys.Ak,1);
Ni = size(sys.Bk,2);
Nd = size(sys.Fk,2);
Nw = size(sys.Gk,2);

sysb.Ns = Ns;
sysb.Ni = Ni;
sysb.Nd = Nd;
sysb.Nw = Nw;

Nsb = Ns*(N+1);
Nib = Ni*(N+1);
Ndb = Nd*(N+1);
Nwb = Nw*(N+1);
Ab = nan(Nsb,Ns);
Bb = nan(Nsb,Nib);
Fb = nan(Nsb,Ndb);
Gb = nan(Nsb,Nwb);
sigk = nan(N+1,1);

Ab(1:Ns,1:Ns) = eye(Ns);
Bb(1:Ns,:) = 0;
Fb(1:Ns,:) = 0;
Gb(1:Ns,:) = 0;

AA = eye(Ns);
AB = Bb(1:Ns,:);
AF = Fb(1:Ns,:);
AG = Gb(1:Ns,:);
for k = 1:N
    
    sigk(k) = sys.Hk(:,:,k)*sys.Swk(:,:,k)*sys.Hk(:,:,k)' + sys.Svk(:,:,k);
    
    rows = (k*Ns+1):(k+1)*Ns;
    colsAB = ((k-1)*Ni+1):k*Ni;
    colsAF = ((k-1)*Nd+1):k*Nd;
    colsAG = ((k-1)*Nw+1):k*Nw;

    Ak = sys.Ak(:,:,k);
    
    % multiply A through input matrices
    for ii = 1:k
        c = (ii-1)*Ni+1:ii*Ni;
        AB(:,c) = Ak*AB(:,c);
    end
    
    for ii = 1:k 
        c = (ii-1)*Nd+1:ii*Nd;
        AF(:,c) = Ak*AF(:,c);
    end
   
    for ii = 1:k
        c = (ii-1)*Nw+1:ii*Nw;
        AG(:,c) = Ak*AG(:,c);
    end
    
    % last columns for inputs/state
    AB(:,colsAB) = sys.Bk(:,:,k);
    AF(:,colsAF) = sys.Fk(:,:,k);
    AG(:,colsAG) = sys.Gk(:,:,k);
    AA = Ak*AA;
    
    Bb(rows,:) = AB;
    Fb(rows,:) = AF;
    Gb(rows,:) = AG;
    Ab(rows,:) = AA;
end
% sigk(N+1) = sys.Hk(:,:,N+1)*sys.Swk(:,:,N+1)*sys.Hk(:,:,N+1)' + sys.Svk(:,:,N+1);

sysb.peakCost = sys.peakCost;
% sysb.alph = sys.alph;
sysb.stateLowerLimits = sys.stateLowerLimits;
sysb.stateUpperLimits = sys.stateUpperLimits;
sysb.controlLowerLimits = sys.controlLowerLimits;
sysb.controlUpperLimits = sys.controlUpperLimits;
sysb.dt = sys.dt;
sysb.sigk = sigk;
sysb.x0 = sys.x0;
sysb.Ab = Ab;
sysb.Bb = Bb;
sysb.Cb = Cb;
sysb.Db = Db;
sysb.Eb = Eb;
sysb.Fb = Fb;
sysb.Gb = Gb;
sysb.Hb = Hb;
sysb.peakCost = sys.peakCost ;
sysb.damageCost = sys.damageCost; 
sysb.energyCost = sys.energyCost; 


