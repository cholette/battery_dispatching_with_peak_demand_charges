function p = normcdf2(x,varargin)

if nargin > 1
    m = varargin{1};
    s = varargin{2};
    x = (x-m)/s;
end
p = 1/2*erfc(-x/sqrt(2));