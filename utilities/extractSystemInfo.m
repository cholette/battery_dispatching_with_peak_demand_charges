function [muk,S] = extractSystemInfo(sys,xk,dk,k,varargin)

muk =  [    sys.Ak(:,:,k)*xk+sys.Fk(:,:,k)*dk 
            sys.Ck(:,:,k)*xk + sys.Ek(:,:,k)*dk];

S = [   sys.Gk(:,:,k)*sys.Swk(:,:,k)*sys.Gk(:,:,k)'      sys.Gk(:,:,k)*sys.Swk(:,:,k)*sys.Hk(:,:,k)'
        sys.Hk(:,:,k)*sys.Swk(:,:,k)*sys.Gk(:,:,k)'      sys.Hk(:,:,k)*sys.Swk(:,:,k)*sys.Hk(:,:,k)' + sys.Svk(:,:,k) ];

% if cond(S) > 1e5
    if nargin > 4
        tol = varargin{1};
    else
        tol = 1e-3;
    end
%    
% end
s = svd(S);
S = S + tol*s(1)*eye(3);