function [xi,phi,s] = setupRBF(lb,ub,nc,type,span,varargin)

% rbf setup
if type == "raw_grid3"
    x = linspace(lb(1),ub(1),nc(1));
    y = linspace(lb(2),ub(2),nc(2));
    z = linspace(lb(3),ub(3),nc(3));
    [X,Y,Z] = ndgrid(x,y,z);
    xi = [X(:) Y(:) Z(:)]';
    
    % based on corners of hypercubes
    dx = diff(x); dx=dx(1);
    dy = diff(y); dy=dy(1);
    dz = diff(z); dz=dz(1);
    s = span*[dx dy dz];
    
elseif lower(type)=="sobol"
    p = sobolset(length(lb));
    nc = prod(nc);
    xi = net(p,nc)';
    xi = xi.*(ub-lb)+lb;
    s = span;
end

if nargin > 5 && ~isempty(varargin{1})% use bias
    phi = @(x) basisFunctions(x,xi,s,"raw",varargin{1},lb,ub);
else % no bias
    phi = @(x) basisFunctions(x,xi,s,"raw",[],lb,ub);
end