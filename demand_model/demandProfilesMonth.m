function traces = demandProfilesMonth(dat,doy1,doy2)
% Extract demand profiles

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
times = dat.timestamps(ind);

% group by year
dt = mode(diff(tod)); % delta t (mode removes uncommon dropped samples and day changes)
nSamples = (doy2-doy1+1)*24/hours(dt);
YR = unique(yr); Ny = length(YR);
traces.demand = nan(nSamples,Ny);
traces.year = YR;
traces.timestamps = NaT(nSamples,Ny);
for ii = 1:Ny
    ind = find( yr==YR(ii));
    traces.demand(:,ii) = demand(ind);
    traces.timestamps(:,ii) = times(ind);
    if ~issorted(traces.timestamps(:,ii) )
        [~,mask] = sort(traces.timestamps(:,ii),'ascend');
        traces.demand(:,ii) = traces.demand(mask,ii);
        traces.timestamps(:,ii) = traces.timestamps(mask,ii);
    end
end