function x = replaceInf(x,M)

ind = find(isinf(x));
x(ind) = sign(x(ind))*M;