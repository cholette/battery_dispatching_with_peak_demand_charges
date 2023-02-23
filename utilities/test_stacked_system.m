clearvars, close all
Ns = 5;
Ni = 5;
Nd = 2;
Nw = 3;
No = 1;
T = 20;
x0 = 100*randn(Ns,1);

sys.Ak = randn(Ns,Ns,T);
sys.Bk = 10*randn(Ns,Ni,T);
sys.Fk = 10*rand(Ns,Nd,T);
sys.Gk = 10*rand(Ns,Nw,T);
sys.Swk = 1*rand(Nw,Nw,T);

sys.Ck = 10*rand(No,Ns,T);
sys.Dk = 10*randn(No,Ni,T);
sys.Ek = 10*randn(No,Nd,T);
sys.Hk = 10*rand(No,Nw,T);
sys.Svk = 1*randn(No,No,T);
sys.peakCost = 0;
sys.damageCost=0;
sys.energyCost=1;
sys.stateLowerLimits=-inf(Ns,1); 
sys.stateUpperLimits=inf(Ns,1);
sys.controlLowerLimits=-inf(Ni,1);
sys.controlUpperLimits=inf(Ni,1);
sys.dt = 1;

U = 100*randn(Ni,T);
W = 100*randn(Nw,T);
V = 100*randn(No,T);

Di = zeros(Nd,T);
for n = 1:Nd
    Di(n,:) = 100*sin( 10*(0:T-1) ) +10*randn(1,T);
end
sys.x0 = x0;

for ii = 1:T
    sys.Swk(:,:,ii) = sys.Swk(:,:,ii)*sys.Swk(:,:,ii)';
    sys.Svk(:,:,ii) = sys.Svk(:,:,ii)*sys.Svk(:,:,ii)';
end

sysb = makeStackedSystem(sys,T-1);

X = 1e7*ones(Ns,T+1);
X(:,1) = x0;
Y = 1e7*ones(No,T);
for t = 1:T
    X(:,t+1) = sys.Ak(:,:,t)*X(:,t) + sys.Bk(:,:,t)*U(:,t) + ...
        sys.Fk(:,:,t)*Di(:,t) + sys.Gk(:,:,t)*W(:,t);
    Y(:,t) = sys.Ck(:,:,t)*X(:,t) + sys.Dk(:,:,t)*U(:,t) +...
        sys.Ek(:,:,t)*Di(:,t) + sys.Hk(:,:,t)*W(:,t) + V(:,t);
end
X(:,end) = [];


Ub = reshape(U,numel(U),1);
Wb = reshape(W,numel(W),1);
Vb = reshape(V,numel(V),1);
Dib = reshape(Di,numel(Di),1);

X_block = sysb.Ab*x0+sysb.Bb*Ub+sysb.Fb*Dib+sysb.Gb*Wb;
Y_block = sysb.Cb*X_block+sysb.Db*Ub+sysb.Eb*Dib+sysb.Hb*Wb+Vb;
X_block = reshape(X_block,[Ns,T]);
Y_block = reshape(Y_block,[No,T]);

figure
tiledlayout(Ns,1)
for n=1:Ns
    plot(X(n,:))
    hold on
    plot(X_block(n,:),'r--')
    nexttile
end
plot(X(n,:))
hold on
plot(X_block(n,:),'r--')

figure
tiledlayout(No,1)
for n=1:No
    plot(Y(n,:))
    hold on
    plot(Y_block(n,:),'r--')
    nexttile
end
plot(Y(n,:))
hold on
plot(Y_block(n,:),'r--')


