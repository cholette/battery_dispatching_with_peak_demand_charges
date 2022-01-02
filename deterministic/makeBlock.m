function XB = makeBlock(varargin)

if ndims(varargin{1}) < 3
    T = varargin{2}; % number of replications
    X = varargin{1};
    XB = cell(1,T);
    [XB{:}]=deal(X);
    XB = blkdiag(XB{:});
elseif ndims(varargin{1})==3 % no need for number of replications
    X = varargin{1};
    T = size(X,3);
    XB = cell(1,T);
    for  ii = 1:T
       XB{ii} = X(:,:,ii); 
    end
    XB = blkdiag(XB{:});
end