function [results,latexTable] = resultsTable(dat)

format = "%.1f $\\pm$ %.1f";
Nadps = size(dat.costADPs,2);
NLookup = size(dat.costLookup,2);
results = strings(Nadps+2,5);
p = dat.peakUC;
cc = dat.costUC;
pc = squeeze(dat.costComponentsUC(:,:,3));
ec = squeeze(dat.costComponentsUC(:,:,1));
dc = squeeze(dat.costComponentsUC(:,:,2));
results(1,1) = sprintf(format,mean(p),std(p));
results(1,2) = sprintf(format,mean(cc),std(cc));
results(1,3) = sprintf(format,mean(pc),std(pc));
results(1,4) = sprintf(format,mean(ec),std(ec));
results(1,5) = sprintf(format,mean(dc),std(dc));

%pk
p = dat.peakPK;
cc = dat.costPK;
pc = squeeze(dat.costComponentsPK(:,:,3));
ec = squeeze(dat.costComponentsPK(:,:,1));
dc = squeeze(dat.costComponentsPK(:,:,2));
results(2,1) = sprintf(format,mean(p),std(p));
results(2,2) = sprintf(format,mean(cc),std(cc));
results(2,3) = sprintf(format,mean(pc),std(pc));
results(2,4) = sprintf(format,mean(ec),std(ec));
results(2,5) = sprintf(format,mean(dc),std(dc));

%mpc
p = dat.peakMPC;
cc = dat.costMPC;
pc = squeeze(dat.costComponentsMPC(:,:,3));
ec = squeeze(dat.costComponentsMPC(:,:,1));
dc = squeeze(dat.costComponentsMPC(:,:,2));
results(3,1) = sprintf(format,mean(p),std(p));
results(3,2) = sprintf(format,mean(cc),std(cc));
results(3,3) = sprintf(format,mean(pc),std(pc));
results(3,4) = sprintf(format,mean(ec),std(ec));
results(3,5) = sprintf(format,mean(dc),std(dc));

%mpc2
p = dat.peakMPC2;
cc = dat.costMPC2;
pc = squeeze(dat.costComponentsMPC2(:,:,3));
ec = squeeze(dat.costComponentsMPC2(:,:,1));
dc = squeeze(dat.costComponentsMPC2(:,:,2));
results(4,1) = sprintf(format,mean(p),std(p));
results(4,2) = sprintf(format,mean(cc),std(cc));
results(4,3) = sprintf(format,mean(pc),std(pc));
results(4,4) = sprintf(format,mean(ec),std(ec));
results(4,5) = sprintf(format,mean(dc),std(dc));

for ii = 1:NLookup
    p = dat.peakLookup(:,ii);
    cc = dat.costLookup(:,ii);
    pc = squeeze(dat.costComponentsLookup(:,ii,3));
    ec = squeeze(dat.costComponentsLookup(:,ii,1));
    dc = squeeze(dat.costComponentsLookup(:,ii,2));
    results(ii+4,1) = sprintf(format,mean(p),std(p));
    results(ii+4,2) = sprintf(format,mean(cc),std(cc));
    results(ii+4,3) = sprintf(format,mean(pc),std(pc));
    results(ii+4,4) = sprintf(format,mean(ec),std(ec));
    results(ii+4,5) = sprintf(format,mean(dc),std(dc));
end

for ii = 1:Nadps
    p = dat.peakADPs(:,ii);
    cc = dat.costADPs(:,ii);
    pc = squeeze(dat.costComponentsADPs(:,ii,3));
    ec = squeeze(dat.costComponentsADPs(:,ii,1));
    dc = squeeze(dat.costComponentsADPs(:,ii,2));
    results(ii+NLookup+4,1) = sprintf(format,mean(p),std(p));
    results(ii+NLookup+4,2) = sprintf(format,mean(cc),std(cc));
    results(ii+NLookup+4,3) = sprintf(format,mean(pc),std(pc));
    results(ii+NLookup+4,4) = sprintf(format,mean(ec),std(ec));
    results(ii+NLookup+4,5) = sprintf(format,mean(dc),std(dc));
end

% create table to export to latex
latexTable = join(results,' & ');
latexTable(:,end) = latexTable(:,end) + "\\";




