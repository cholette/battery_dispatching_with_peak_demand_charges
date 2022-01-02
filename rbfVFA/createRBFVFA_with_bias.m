function rbfVFA = createRBFVFA_with_bias(lb,ub,nc,T,type,span)
% Initialize a radial basis function value function approximaiton

rbfVFA.lb = lb;
rbfVFA.ub = ub;
disp(['Initializing RBFVFA from scratch for ',num2str(length(lb)) ' dimensional state with ',...
    num2str(prod(nc)),' basis functions (from zero)'])
[rbfVFA.xi,rbfVFA.phi,rbfVFA.sigma] = setupRBF(lb,ub,nc,type,span,'bias');
nc = max(size(rbfVFA.xi));
rbfVFA.theta = zeros(nc+1,T+1);
rbfVFA.Brls = [];
rbfVFA.nc = nc;
rbfVFA.type = type;
rbfVFA.span = span;
rbfVFA.T = T;

