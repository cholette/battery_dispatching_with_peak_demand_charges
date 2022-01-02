function traces = demandProfiles(dat,doy1,doy2,varargin)
% Extract daily demand profiles between two days of the 
% year.

% find day of Year for each entry
doy = day(dat.timestamps,'dayofyear');
tod = timeofday(dat.timestamps);
yr = year(dat.timestamps);
traces.dt = mode(diff(dat.timestamps));

% adjust DOY for leapyear
doy = doy - (leapyear(yr)&doy>60);

% extract samples within date range
ind = find(doy>=doy1 & doy<=doy2);
doy = doy(ind);
tod = tod(ind);
yr = yr(ind);
demand = dat.kW(ind);

if nargin>3 && strcmpi(varargin{1},'flat')
    traces.doy=doy;
    traces.tod=tod;
    traces.year = yr;
    traces.demand = demand;
    return
end

% get number of samples in a full day
dt = mode(diff(tod)); % delta t (mode removes uncommon dropped samples and day changes)
nSamples = duration(24,0,0)/dt;


% extract daily traces
YR = unique(yr); Ny = length(YR);
DD = unique(doy); Nd = length(DD);
traces.demand = nan(nSamples,Nd,Ny);
traces.year = YR;
traces.tod = min(tod):dt:max(tod);
for ii = 1:Ny
   for jj = 1:Nd
       ind = find( yr==YR(ii) & doy==DD(jj) ); 
       if ~isempty(ind) && length(ind)==nSamples % only complete days
           traces.demand(:,jj,ii) = demand(ind);
       end
   end
end

if nargin>3 && strcmpi(varargin{1},'group')
   Ntod = length(traces.tod);
   D = traces.demand;
   D = reshape(D,[Ntod,Nd*Ny]);
   traces.demand = D;
   traces.doy = repmat(DD',[1,Ny]);
end

