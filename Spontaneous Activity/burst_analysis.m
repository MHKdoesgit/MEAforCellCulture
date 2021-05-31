
function out = burst_analysis(spk, recduration)

plotting    =   false;

ncells      =   numel(spk);
binlength   =   100/1e3;
recdur      =   ceil(recduration/5)*5;

thresh.rasterburstduration      =       2;
thresh.rasterburstend           =       1;
thresh.psthburstonset           =       0.15;
thresh.psthburstduration        =       0.5;
thresh.psthpeakratio            =       5;
thresh.meanspkperbursts         =       5;
thresh.maxrateratio             =       1.5;

% get all the spikes in a matrix
spkall = CelltoMatUE(spk);


[raster_burststart,raster_burstend, raster_spkperburst] = deal(cell(ncells,1));
[diffspk, burstidx] = deal(cell(ncells,1));

for jj = 1: ncells
    s = spk{jj};
    diffspk{jj} = diff([0;s]);
    
    burstidx{jj} = diffspk{jj} > thresh.rasterburstduration;
    raster_burststart{jj} = s(burstidx{jj});
    raster_burstend{jj} = s(burstidx{jj}) + thresh.rasterburstend;
    
    spkpb = zeros(sum(burstidx{jj}), 1);
    for kk = 1: sum(burstidx{jj})
        spkpb(kk) = sum(s > raster_burststart{jj}(kk) & s <= raster_burstend{jj}(kk));
    end
    raster_spkperburst{jj} = spkpb;
end


% three thresholds are defined to select the bursty cells.
% 1- more than 1 spikes in any detected burst
% 2- At least 5 spikes across all the detected bursts
% 3- less than 1.5 Hz firing rate to avoid tonic cells with not bursts
burstycells = all([not(cellfun('isempty',raster_spkperburst)) ,...
    cellfun('length',raster_spkperburst) > thresh.meanspkperbursts,...
    cellfun('length',spk)/recdur < thresh.maxrateratio],2);

spkburstycells = spkall(burstycells,:);

% Now we construct a psth from bursty cells
binvec = linspace(0, recdur,(recdur/binlength));
% psth = histcounts(spkburstycells, binvec);        % histcounts suck!
psth = histc(spkburstycells',binvec) * (1/binlength);   %#ok
% for all the cells
psthall = histc(spkall',binvec) * (1/binlength);   %#ok

% calculate reliablitiy Rsq for bursty PSTH and all PSTHs
psthodds = mean( psth(:,1:2:end), 2); 
psthevens = mean( psth(:,2:2:end), 2);
[rsq, pearsoncoeff, explainedvar] = calcRsqPearsonCoeff(psthodds,psthevens);

psthodds = mean( psthall(:,1:2:end), 2); 
psthevens = mean( psthall(:,2:2:end), 2);
[rsq_all, pearsoncoeff_all, explainedvar_all] = calcRsqPearsonCoeff(psthodds,psthevens);

psth = mean (psth,2) ;
psthall = mean (psthall,2) ;

% this is to take the average of the three biggest peaks to set the threshold
sortedpsth = sort(psth,'descend');
peakthresh = mean(sortedpsth(1:3))/thresh.psthpeakratio;

% Now find the peaks that are bigger than the threshold
[burstpeaks, burstlocs] = findpeaks(psth,'MinPeakHeight',peakthresh,'MinPeakDistance',2*thresh.psthburstduration/binlength);
burstlocstvec = binvec(burstlocs)';

burststart = burstlocstvec - thresh.psthburstonset;
burstend = burstlocstvec + thresh.psthburstduration;
numbursts = numel(burstpeaks);

bps = cell(numbursts,1);
[numspkperburst, burstduration] = deal(zeros(ncells, numbursts));
burstspks = cell(ncells,1);

for jj = 1:ncells
    for kk = 1: numbursts
        bps{kk} = spk{jj}(spk{jj} > burststart(kk) & spk{jj} <= burstend(kk));
        numspkperburst(jj,kk) = numel(bps{kk});
    end
    burstspks{jj} = CelltoMatUE(bps)';
    if not(isempty(burstspks{jj}))
    burstduration(jj,:) = (nanmax(burstspks{jj}) - nanmin(burstspks{jj}))*1000;
    end
end
numspkperburst = numspkperburst(burstycells,:);
burstduration = burstduration(burstycells,:);
burstduration(burstduration==0) = NaN;
burstspks = burstspks(burstycells,:);


%    burstycells = all([not(cellfun('isempty',spkperburst)) , (sum(numspkperburst,2) > thresh.meanspkperbursts ),...
%        cellfun('length',spk)/recdur < thresh.maxrateratio],2);

indices.firingRate = (cellfun('length',spk(burstycells))/recduration);
indices.burstpermin = numbursts/(recduration/60);
indices.ibi = diff(burstlocstvec);
indices.intraburstfrq =  nanmean((numspkperburst./burstduration)*1e3,2); % back to Hz here
indices.spkinburstpercent = sum(numspkperburst,2)./cellfun('length',spk(burstycells));
indices.burstduration = nanmean(burstduration,2);
indices.Rsq = rsq;
indices.pearsoncoeff = pearsoncoeff;
indices.explainedvar = explainedvar;
indices.Rsq_allcells = rsq_all;
indices.pearsoncoeff_allcells = pearsoncoeff_all;
indices.explainedvar_allcells = explainedvar_all;
indices.totalspikes = cellfun('length',spk);
indices.totalburstyspikes = cellfun('length',spk(burstycells));
indices.allchanFiringRate = cellfun('length',spk)/recduration;
indices.allmeanISI = cellfun(@mean,cellfun(@diff,spk,'un',0));
indices.meanISIbursty = cellfun(@mean,cellfun(@diff,spk(burstycells),'un',0));


mx = max(burstduration,[],'all');
if mx > thresh.psthburstduration
    try
        mx = max(burstduration(~isoutlier(burstduration,'gesd')),[],'all');
    catch
        mx = max(burstduration(~isoutlier(burstduration)),[],'all');
    end
end
raspsthallbursts = raspsthperburst(burststart, spk(burstycells), mx);


out.rastersall = spkall;
out.rastersbursty = spkburstycells;
out.psthall = psthall;
out.psth = psth;
out.psthtime = binvec;
out.burstycellindex = burstycells;
out.raspsthallbursts = raspsthallbursts;
out.peakthreshold = peakthresh;
out.burstpeaks = burstpeaks;
out.burstlocs = burstlocs;
out.burstpeaktime = burstlocstvec;
out.burststart = burststart;
out.burstend = burstend;
out.numbursts = numbursts;
out.numspksperburst = numspkperburst;
out.burstduration = burstduration;
out.burstrasters = burstspks;
out.indices = indices;
out.thresholds = thresh;
out.binlength = binlength;
out.recordingdur = recduration;

%%
if plotting
    figure('pos',[1 50 1900 950],'color','w');
    cellid = 29 ;
    yAx = size(spkburstycells,1)+3;
    subplot_tight(2,1,1,0.04)
    patch([burststart, burstend, burstend,burststart]',repmat([0 0 yAx, yAx],size(burststart,1),1)',...
        'y','edgecolor','none','facealpha',0.75);
    hold on;
    rasterPlotter(spkburstycells',[],'k');   hold on;
    %plot([burststart burststart],[0 yAx],'r--');
    %plot([burstend burstend],[0 yAx],'g--');
    rasterPlotter(spkburstycells(cellid,:)',[],'m',cellid);
    title(['cellid: ', num2str(cellid)]);
    axis([0 binvec(end-1) 0 yAx]);
    set(gca,'ticklength',[0.0025 0.0025]);
    
    
    subplot_tight(2,1,2,0.04)
    yAx = ceil(max(psth)/5)*5;
    plot(binvec(1:end-1), psth,'k-',burstlocstvec,burstpeaks,'ro');    %pbaspect([3,1,1]);
    hold on;
    plot([burststart burststart],[0 yAx],'g--');
    %plot([bend bend],[0 size(spkburstycells,1)+3],'r--');
    axis([0 binvec(end-1) 0 yAx]);
    set(gca,'ticklength',[0.0025 0.0025]);        box off;
    
end

%%
%        %[~,score] = pca(d,'NumComponents',2);
%
%        [y,x] = hist(d,1000);
%        G = spline(x,y);
%
%        % Find optimum curve fit
%        P0 = [% mu  S    A
%            -0.05  0.1   1;  % (some rough initial estimate)
%            0.25  0.5  0.5];
%        P = fminunc(@(P) Obj(P, x,G), P0); % refine the estimate
%
%
%        % REPRODUCED DATA
%        %P(:,1:2).'
%
%        figure, clf, hold on
%        plot(x, P(1,3)*Gaussian(P(1,1),P(1,2),x) + P(2,3)*Gaussian(P(2,1),P(2,2),x))
%        plot(x, ppval(G,x),'r.', 'MarkerSize', 1)
%
%
%        km = kmeans(d,2);
%
%        GMModel = fitgmdist(d,2);
%
%



end


function out = raspsthperburst(ft, spks, maxbrdur)

% stim = ft(1:steps:numel(ft));
binLength = 10/1e3;
stimround = ceil(maxbrdur/1e3); % round up to closest second
if stimround == 0, stimround = 0.1; end % to avoid crashing by zero
binVec = linspace(0,stimround,stimround/binLength);

npulses = numel(ft);
ncells = numel(spks);

rasters = cell(npulses,1);
psth = zeros(npulses, length(binVec));

% ft= [ft(1:end);ft(end)+stimround]; 


for ii = 1: npulses
    ras = cell(1,ncells);
    for jj = 1: ncells
        ras{jj} = spks{jj}(spks{jj} > ft(ii) & spks{jj} <= ft(ii)+stimround) - ft(ii);
    end
    rasters{ii} = CelltoMatUE(ras);
    
    if size(rasters{ii},2) < 2
        try
            psth(ii,:) = (histc(rasters{ii}',binVec)/size(rasters{ii},1))*(1/binLength);  %#ok to fix situation with 1 spike
        catch
            psth(ii,:) = mean(histc(rasters{ii}',binVec),2)*(1/binLength);   %#ok
        end
    else
        psth(ii,:) = mean(histc(rasters{ii}',binVec),2)*(1/binLength);   %#ok must transpose contspk for histc to work properly.
    end
end
psthall = mean(psth,1)';


[lat, pk] = psthPeakLatency(psthall,[5 300],5, 1000);

[psthdctime, dcval,dctvec] = psthDecayTime(psthall, lat, 1000);

psthodds = mean( psth(1:2:end,:), 1); 
psthevens = mean( psth(2:2:end,:), 1);

[rsq, pearsoncoeff, explainedvar] = calcRsqPearsonCoeff(psthodds,psthevens);


out.rasters = rasters;
out.psth = psthall;
out.psthperburst = psth;
out.psthtime = binVec;
out.peak = pk;
out.latency = lat;
out.psthdecayindex = psthdctime;
out.decayval = dcval;
out.decaytvec = dctvec;
out.rsq = rsq;
out.pearsoncoeff = pearsoncoeff;
out.explainedvar = explainedvar;
out.binlength = binLength;
out.stimdur = stimround;


end




function val = Obj(P, x,F)

G = zeros(size(x));
for ii = 1:size(P,1)
    
    mu = P(ii,1);    % mean
    sigma = P(ii,2); % std. deviation
    A = P(ii,3);     % "amplitude"
    
    G = G + A/sigma/sqrt(2*pi) * exp(-(x-mu).^2/2/sigma^2);
    
end

val = sum((G-ppval(F,x)).^2);

end

% just a function for plotting
function G = Gaussian(mu,sigma,x)
G = 1/sigma/sqrt(2*pi) * exp(-(x-mu).^2/2/sigma^2);
end




function [gaussPdfi, binneddata, xaxis, gmdist] = fitgaussianmixture(nl, Numgaussian, nbins)

%xaxis = transpose(min(nl):range(nl)/(nbins-1):max(nl));
[binneddata, xaxis]  = hist(nl,nbins);
xaxis = transpose(xaxis);
% pre-allocation
mu     = zeros(Numgaussian,1);
sigma  = zeros(Numgaussian,1);
weight = zeros(Numgaussian,1);
gaussPdfi = zeros(nbins,Numgaussian);
%gaussPdf  = zeros(nbins,1);

% fit
options = statset('Display','final');
obj = gmdistribution.fit(nl,Numgaussian,'Options',options);
gaussPdf = pdf(obj,xaxis);
A = sum(gaussPdf);
gaussPdf = gaussPdf/A;

% separating N Gaussians
for n = 1:Numgaussian
    mu(n)          = obj.mu(n);
    sigma(n)       = sqrt(obj.Sigma(1,1,n));
    weight(n)      = obj.PComponents(n);
    gaussPdfi(:,n) = weight(n)*normpdf(xaxis,mu(n),sigma(n))/A;
end

% yrangediff = (max(binneddata)/max(gaussPdfi(:)));
% gaussPdfi = gaussPdfi * yrangediff;
gmdist.gaussPdf = gaussPdf;
gmdist.gaussPdfi = gaussPdfi;
gmdist.mu = mu;
gmdist.sigma = sigma;
gmdist.weight = weight;
gmdist.obj = obj;



plot(xaxis,gaussPdfi,'k','linewidth',1);


end


