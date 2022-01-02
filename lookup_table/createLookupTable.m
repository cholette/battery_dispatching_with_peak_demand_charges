function lookupVFA = createLookupTable(N,LB,UB,T)

nd = length(LB);
if nd ~= length(N)
    error('Number of cells of N must be equal to the number of dimensions (length(LB) & length(UB))')
end

% create uniformly-spaced grid in each dimension
XI = cell(1,nd);
lookupVFA.dx = nan(1,nd);
for ii = 1:nd
    XI{ii} = linspace(LB(ii),UB(ii),N(ii));
    lookupVFA.dx(ii) = (UB(ii)-LB(ii))/(N(ii)-1);
end
[XI{:}] = ndgrid(XI{:});

Npts = numel(XI{1});
lookupVFA.xi = zeros(Npts,nd);
for ii = 1:nd
    lookupVFA.xi(:,ii) = XI{ii}(:);
end

lookupVFA.values = zeros(Npts,T);
lookupVFA.Ng = Npts;


