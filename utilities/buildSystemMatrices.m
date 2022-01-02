function [sys,dk] = buildSystemMatrices(base,demandModel,t)
% builds system matrices from excel file and the demand PAR(1) model

% load in matrices
A = base.A;
B = base.B;
C = base.C;
D = base.D;
F = base.F;
E = base.E;
G = base.G;
H = base.H;

N = length(t);
Ns = size(A,1);
Ni = size(B,2);
No = size(C,1);
Nd = size(F,2);
Nw = size(G,2);
sys.Ak = nan(Ns,Ns,N);
sys.Bk = nan(Ns,Ni,N);
sys.Ck = nan(No,Ns,N);
sys.Dk = nan(No,Ni,N);
sys.Ek = nan(No,Nd,N);
sys.Fk = nan(Ns,Nd,N);
sys.Gk = nan(Ns,Nw,N);
sys.Hk = nan(No,Nw,N);
sys.Swk = nan(Nw,Nw,N);
sys.Svk = nan(No,No,N);
dk = zeros(1,N);

for nn = 1:N
    % get Ak
    sys.Ak(:,:,nn) = A;
    cc = find(demandModel.tod==timeofday(t(nn)));
    sys.Ak(1,1,nn) = demandModel.ar(1,cc);
    
    % get other system matrices, most of which are
    % constant.
    sys.Bk(:,:,nn) = B;
    sys.Ck(:,:,nn) = sys.Ak(1,:,nn);
    sys.Dk(:,:,nn) = D;
    sys.Ek(:,:,nn) = E;
    sys.Fk(:,:,nn) = F;
    sys.Gk(:,:,nn) = G;
    sys.Hk(:,:,nn) = H;
    
    % disturbance characteristics
    sys.Swk(:,:,nn) = demandModel.std(cc)^2;
    sys.Svk(:,:,nn) = zeros( size(sys.Ck,1),1);
    dk(nn) = demandModel.ar(2,cc);
end
