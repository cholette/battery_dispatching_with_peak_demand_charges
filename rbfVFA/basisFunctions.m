function PHI = basisFunctions(x,xi,s,type,varargin)

nc = size(xi,2);
if nargin < 5 || isempty(varargin{1})
    PHI = nan(nc,1);
elseif varargin{1}=="bias"
    PHI = nan(nc+1,1);
    PHI(end) = 1;
else
    error('bias option not recognized')
end

if isscalar(s)
    n = size(xi,1);
    s = s*ones(1,n);
end

% d = (repmat(x,[1,nc])-xi);
% w2 = exp( -sum(d.^2)/2/s(1)^2);
d2 = pdist2(xi',x','seuclidean',s);
w = exp( -1/2*d2.^2);
switch type
    case "normalized"
        PHI(1:nc) = w./sum(w);
    case "raw"
        PHI(1:nc) = w;
end

