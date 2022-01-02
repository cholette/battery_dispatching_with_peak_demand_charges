function compareCosts(dat,monthStr,modelInds,nbins,fontSize,ADPNames)

% COSTS
nModels = length(modelInds);
hf = figure("Name",monthStr+" Costs");
subplot(nModels+2,1,1);
histogram(dat.costUC,nbins,'normalization','probability'); title('Uncontrolled')
subplot(nModels+2,1,2); histogram(dat.costPK,nbins,'normalization','probability'); title('Perfect')
MEANS = [mean(dat.costUC),mean(dat.costPK)];%,mean(datTrainingSim.costADP),mean(dat_search.costH)];
STDS = [std(dat.costUC),std(dat.costPK)];%,std(datTrainingSim.costADP),std(datTrainingSim.costH)];
for ii = 1:nModels
    COSTS = dat.costADPs(:,modelInds(ii));
    subplot(2+nModels,1,2+ii); histogram(COSTS,nbins,'normalization','probability'); 
    title(ADPNames(ii))
    MEANS(end+1) = mean(COSTS);
    STDS(end+1) = std(COSTS);
end
hf = formatPlots(hf,fontSize,"","Probability");
xlabel(hf.Children(1),"Cost (AUD)")
hf.Units = "normalized";
hf.Position = [0.1 0.2 0.25 0.6];
meansOnPlots(hf,fliplr(MEANS),fliplr(STDS)); % order of children is front to back

% Peaks
hf = figure("Name",monthStr+" Peaks");
subplot(nModels+2,1,1); 
histogram(dat.peakUC,nbins,'normalization','probability'); title('Uncontrolled')
subplot(nModels+2,1,2); 
histogram(dat.peakPK,nbins,'normalization','probability'); title('Perfect')
MEANS = [mean(dat.peakUC),mean(dat.peakPK)];
STDS = [std(dat.peakUC),std(dat.peakPK)];
for ii = 1:nModels
    PEAKS = dat.peakADPs(:,modelInds(ii));
    subplot(2+nModels,1,2+ii); histogram(PEAKS,nbins,'normalization','probability'); 
    title(ADPNames(ii))
    MEANS(end+1) = mean(PEAKS);
    STDS(end+1) = std(PEAKS);
end
hf = formatPlots(hf,fontSize,"","Probability");
xlabel(hf.Children(1),"Peak Demand (kWh)")
hf.Units = "normalized";
hf.Position = [0.4 0.2 0.25 0.6];
meansOnPlots(hf,fliplr(MEANS),fliplr(STDS)); % order of children is front to back

% Energy Costs
ECPK = squeeze(dat.costComponentsPK(:,:,1));
ECUC = squeeze(dat.costComponentsUC(:,:,1));
ECADPs = squeeze(dat.costComponentsADPs(:,:,1));

hf = figure("Name",monthStr+" Energy Costs");
subplot(nModels+2,1,1); 
histogram(ECUC,nbins,'normalization','probability'); title('Uncontrolled')
subplot(nModels+2,1,2); 
histogram(ECPK,nbins,'normalization','probability'); title('Perfect')
MEANS = [mean(ECUC),mean(ECPK)];
STDS = [std(ECUC),std(ECPK)];
for ii = 1:nModels
    ECA = ECADPs(:,modelInds(ii));
    subplot(2+nModels,1,2+ii); histogram(ECA,nbins,'normalization','probability'); 
    title(ADPNames(ii))
    MEANS(end+1) = mean(ECA);
    STDS(end+1) = std(ECA);
end
hf = formatPlots(hf,fontSize,"","Probability");
xlabel(hf.Children(1),"Energy Cost (AUD)")
hf.Units = "normalized";
hf.Position = [0.4 0.2 0.25 0.6];
meansOnPlots(hf,fliplr(MEANS),fliplr(STDS)); % order of children is front to back

% Peak Costs
PCPK = squeeze(dat.costComponentsPK(:,:,3));
PCUC = squeeze(dat.costComponentsUC(:,:,3));
PCADPs = squeeze(dat.costComponentsADPs(:,:,3));

hf = figure("Name",monthStr+" Peak Costs");
subplot(nModels+2,1,1); 
histogram(PCUC,nbins,'normalization','probability'); title('Uncontrolled')
subplot(nModels+2,1,2); 
histogram(PCPK,nbins,'normalization','probability'); title('Perfect')
MEANS = [mean(PCUC),mean(PCPK)];
STDS = [std(PCUC),std(PCPK)];
for ii = 1:nModels
    PCA = PCADPs(:,modelInds(ii));
    subplot(2+nModels,1,2+ii); histogram(PCA,nbins,'normalization','probability'); 
    title(ADPNames(ii))
    MEANS(end+1) = mean(PCA);
    STDS(end+1) = std(PCA);
end
hf = formatPlots(hf,fontSize,"","Probability");
xlabel(hf.Children(1),"Peak Cost (AUD)")
hf.Units = "normalized";
hf.Position = [0.4 0.2 0.25 0.6];
meansOnPlots(hf,fliplr(MEANS),fliplr(STDS)); % order of children is front to back

% Damage Costs
DCPK = squeeze(dat.costComponentsPK(:,:,2));
DCUC = squeeze(dat.costComponentsUC(:,:,2));
DCADPs = squeeze(dat.costComponentsADPs(:,:,2));

hf = figure("Name",monthStr+" Damage Costs");
subplot(nModels+2,1,1); 
histogram(DCUC,nbins,'normalization','probability'); title('Uncontrolled')
subplot(nModels+2,1,2); 
histogram(DCPK,nbins,'normalization','probability'); title('Perfect')
MEANS = [mean(DCUC),mean(DCPK)];
STDS = [std(DCUC),std(DCPK)];
for ii = 1:nModels
    DCA = DCADPs(:,modelInds(ii));
    subplot(2+nModels,1,2+ii); histogram(DCA,nbins,'normalization','probability'); 
    title(ADPNames(ii))
    MEANS(end+1) = mean(DCA);
    STDS(end+1) = std(DCA);
end
hf = formatPlots(hf,fontSize,"","Probability");
xlabel(hf.Children(1),"Damage Cost (AUD)")
hf.Units = "normalized";
hf.Position = [0.4 0.2 0.25 0.6];
meansOnPlots(hf,fliplr(MEANS),fliplr(STDS)); % order of children is front to back