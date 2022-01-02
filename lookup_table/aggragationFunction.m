function [xg,indg,N] = aggragationFunction(X,lookup,varargin)

dx = lookup.dx;
xi = lookup.xi;
N = size(X,1);
if nargin<3
    method = "pdist2";
else
    method = varargin{1};
end
switch method
    case "pdist2"
        [~,indg] = pdist2(xi,X,'euclidean','Smallest',1);
        xg = xi(indg,:);
    case "loop"
        xg = nan(size(X));
        indg = nan(N,1);
        for ii = 1:size(X,1)
            x = X(ii,:);
            [~,indg(ii)] = min( sum((x-xi).^2,2) );
            xg(ii,:) = xi(indg(ii),:);
        end
    case "edges"
        lb = xi-dx/2;
        ub = xi+dx/2;
        lb(lb<min(xi)) = -inf;
        ub(ub>max(xi)) = inf;   
        X2 = zeros(1,size(X,2),size(X,1));
        X2(1,:,:) = X';
        temp = all(X2<=ub & X2>=lb,2);
        [indg,~] = find(squeeze(temp));
        xg = xi(indg,:);
end
e = (0:length(xi))+0.5;
e(1) = -inf; e(end) = inf;
N = histcounts(indg,e);
