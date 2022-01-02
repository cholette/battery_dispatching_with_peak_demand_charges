function [rbfVFA,phisOut] = estimateRBFParameters(k,X,rbfVFA,Vk,param,plotFlag)

% pre-compute activations
disp("Computing activations & estimating coefficients ...")
Ns = length(X);
if iscell(rbfVFA.theta)
    ncoeff = size(rbfVFA.theta{k},1);
    nc = rbfVFA.nc(k);
else
    ncoeff = size(rbfVFA.theta,1);
    nc = rbfVFA.nc;
end
phisOut = zeros(Ns,ncoeff);
theta = zeros(ncoeff,1);

for ii = 1:Ns
    if iscell(rbfVFA.phi)
        phisOut(ii,:) = rbfVFA.phi{k}(X(:,ii))';
    else
        phisOut(ii,:) = rbfVFA.phi(X(:,ii))';
    end
end

% handle bias
if ncoeff > nc % bias
    b = mean(Vk);
    theta(end) = b;
    Vk = Vk - b;
    ncoeff = ncoeff - 1;
    phis = phisOut(:,1:nc);
else
    phis = phisOut;
end
   
% remove low activations
tolActiv = param(2);
ind = find(max(phis)>tolActiv);
phis2 = phis(:,ind);
disp(num2str(length(ind))+" centers activated above "+num2str(tolActiv))

% setup regularization params
tol = 16*eps;
lams = logspace(log10(tol),0,100);

c = cvpartition(Ns,'HoldOut',param(1));
trmask = training(c,1);
temask = test(c,1);
[U,s,V] = csvd(phis2(trmask,:));
lams = lams*s(1);
x_lam = tikhonov(U,s,V,Vk(trmask),lams);
err = repmat(Vk(temask),1,length(lams)) - phis2(temask,:)*x_lam;
mse = sum(err.^2);
[~,ropt] = min(mse);

% last fitting
[U,s,V] = csvd(phis2);
[temp,rho,eta] = tikhonov(U,s,V,Vk,lams(ropt));
theta(ind) = temp;

if plotFlag
    fh = findobj( 'Type', 'Figure', 'Name', "L-curve" );
    if isempty(fh)
        figure("Name","L-curve")
    else
        figure(fh)
    end
    l_curve(U,s,Vk,"Tikh",plotFlag)
    hold on
    loglog(rho,eta,'r*')
    hold off
end


if iscell(rbfVFA.theta)
    rbfVFA.theta{k} = theta;
else
    rbfVFA.theta(:,k) = theta;
end