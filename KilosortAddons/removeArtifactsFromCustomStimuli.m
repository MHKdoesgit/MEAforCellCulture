

function plsout = removeArtifactsFromCustomStimuli(filetowrite, h5dat, chanMap, nsamples, maxSamples, savingpath, varargin)

if nargin > 6, stimpara = varargin{1}; else, stimpara.stimulus = NaN; end

fprintf('extracting electrical pulses from the raw data...');
% sort out the saving path first
sn = extractBefore(savingpath,'.mat');
% chalabels = cellfun(@str2double, dat.Info.Label);
% stimlabels = cellfun(@str2double, eventstream.Info.SourceChannelIDs);
% find(ismember(chalabels, stimlabels))

ops.fs = h5dat.getSamplingRate;
ops.spkTh           = -50;     % spike threshold in standard deviations (4), -125 for old method
ops.loc_range       = [25 0];  % ranges to detect peaks; plus/minus in time and channel ([3 1])
ops.long_range      = [60  0]; % ranges to detect isolated peaks ([30 6])
ops.nt0             = 351;
ops.nt0min          = 151;
ops.goodchnum       = 3;
ops.minexpectedpulse = 5;


Nchunk              =    ceil(nsamples/maxSamples);
[pulsetemplate, pulsestartstop, pulselocindex, besttemploc]     =   deal(cell(1,Nchunk));
spkth               =    zeros(Nchunk,1);

for iChunk  =   1 : Nchunk
    offset          =   max(0, (maxSamples * (iChunk-1)));
    sampstoload     =   min(nsamples-offset,maxSamples);
    lastidxtoload   =   offset+sampstoload; % added by MHK to avoid end index crashes
    if lastidxtoload > size(h5dat.ChannelDataTimeStamps,2)
        lastidxtoload   =    size(h5dat.ChannelDataTimeStamps,2);
    end
    cfg.window  =   double([h5dat.ChannelDataTimeStamps(offset+1) h5dat.ChannelDataTimeStamps(lastidxtoload)]) / 1e6;
    % to convert from microseconds to sec for more info check McsHDF5.TickToSec
    
    dat         =   h5dat.readPartialChannelData(cfg);
    % rearrange channel to the channel map
    chdat       =    dat.ChannelData(chanMap + 1,:);
    chtvec      =    dat.ChannelDataTimeStamps;
    [pulsetemplate{iChunk}, pulsestartstop{iChunk}, pidx, tidx, spkth(iChunk)] = get_electrical_pulses(chdat, chtvec, ops);
    pulselocindex{iChunk}   =   offset + pidx;
    besttemploc{iChunk}     =   offset + tidx;
    
end

chunkswithpulse = not(cellfun(@(x)all(isnan(x)),besttemploc)); % to deal with chunk that have no pulse

[npulsese, ~,nChan]     =   cellfun(@size,pulsetemplate(chunkswithpulse));
nChan                   =   unique(nChan);
npsum                   =   sum(npulsese);
allpulsetemplate        =   zeros(npsum, ops.nt0, nChan);
[allpulsestartstop, allpulselocindex] = deal(zeros(npsum, ops.nt0));
allbesttemploc          =   zeros(npsum,1);
np                      =   cumsum([1,npulsese]);

for iChan = 1:sum(chunkswithpulse)
    nporder = np(iChan):np(iChan+1)-1;
    allpulsetemplate(nporder,:,:)   =   pulsetemplate{iChan};
    allpulsestartstop(nporder,:)    =   pulsestartstop{iChan};
    allpulselocindex(nporder,:)     =   pulselocindex{iChan};
    allbesttemploc(nporder)         =   besttemploc{iChan};
end

% first get electrical image by averaging all the pulses for each channel
meantemplates   =   squeeze(mean(allpulsetemplate,1));
% next claculate the eucleadian norm of the mean templates
meantempnorm    =   meantemplates ./ sqrt(sum(meantemplates.^2,1));
% now get the norm of each individual pulse
alltempnorm     =  (sqrt(sum(allpulsetemplate.^2,2)));
% finally scale each average pulse to norm of the individual pulses
meantosubtract  =   alltempnorm .* reshape(meantempnorm,1,ops.nt0, nChan);

% figure
% plot(squeeze(allpulsetemplate(1,:,1)))
% hold on
% plot(squeeze(meantosubtract(1,:,1)));
% plot(squeeze(allpulsetemplate(1,:,1))-squeeze(meantosubtract(1,:,1)))
% plot(meantemplates(:,1));


fprintf('%d pulses found, fuck yeah...\n',sum(npulsese));

% this is to test for all the channels at the same time

%  chdatall = single (h5dat.ChannelData);
%  chdatall    =    chdatall(chanMap + 1,:);
%  chtvecall = h5dat.ChannelDataTimeStamps;
% 
%  for ii = 1: nChan
%      thischdat = chdatall(ii, :);
%      m = single (squeeze(meantosubtract(:,:,ii)));
%      thischdat(allpulselocindex) =  thischdat(allpulselocindex) - m ;
%      chdatall(ii,:) = thischdat;
% 
%  end
% chdatallraw = h5dat.ChannelData;
% chdatallraw    =    chdatallraw(chanMap + 1,:);
% %
%  for ii = 1:60
%      subplot(2,1,1);
%      plot(chtvecall,chdatallraw(ii,:)); title(ii); axis tight; ax = gca;
%      subplot(2,1,2);
%      plot(chtvecall,chdatall(ii,:),'r'); set(gca,'YLim',ax.YLim, 'XLim',ax.XLim);
%      pause;
%  end


%  refchanidx = find(chanMap == find(contains(h5dat.Info.Label,'Ref')))+1;
%
if Nchunk > 10 % for more than 10 chunk randomly plot 10 of the chunks
    pltdatflag = false(1,Nchunk);
    pltdatflag(randperm(Nchunk,10)) = true;
else
    pltdatflag = true(1,Nchunk);
end
% tsfolder = extractBefore(savingpath,'electrical_stimuli');


for iChunk  =   1 : Nchunk
    offset          =   max(0, (maxSamples * (iChunk-1)));
    sampstoload     =   min(nsamples-offset,maxSamples);
    lastidxtoload   =   offset+sampstoload; % added by MHK to avoid end index crashes
    if lastidxtoload > size(h5dat.ChannelDataTimeStamps,2)
        lastidxtoload   =    size(h5dat.ChannelDataTimeStamps,2);
    end
    cfg.window  =   double([h5dat.ChannelDataTimeStamps(offset+1) h5dat.ChannelDataTimeStamps(lastidxtoload)]) / 1e6;
    % to convert from microseconds to sec for more info check McsHDF5.TickToSec
    
    dat         =   h5dat.readPartialChannelData(cfg);
    chdat       =   single(dat.ChannelData(chanMap + 1,:));
    chtvec      =    dat.ChannelDataTimeStamps;
    chunkpulseidx = pulselocindex{iChunk} - offset;
    
    
    for iChan = 1: nChan
        if isnan(chunkpulseidx), continue; end
        thischdat   =   chdat(iChan, :);
        mntempl   =   single (squeeze(meantosubtract(np(iChunk):np(iChunk+1)-1,:,iChan)));
        thischdat(chunkpulseidx)    =   thischdat(chunkpulseidx) - mntempl ;
        
        % first plotting before over writing data
        if pltdatflag(iChunk)
            plstempl = single(pulsetemplate{iChunk}(:,:,iChan));
            plotelectstim(chdat(iChan,:), thischdat, chtvec, mntempl, meantemplates(:,iChan),  plstempl ,ops, [iChan,iChunk], sn);
        end
        
        chdat(iChan,:)  =   thischdat;
        
    end
    
    chdattowrite = int16(chdat);
    fwrite(filetowrite, chdattowrite, 'int16');
    
end

plsout.pulsetemplates       =       allpulsetemplate;
plsout.pulsestartstop       =       allpulsestartstop;
plsout.pulseindices         =       allpulselocindex;
plsout.pulsedetectedat      =       allbesttemploc;
plsout.meantemplates        =       meantemplates;
plsout.mean2subtract        =       meantosubtract;
ops.numpulses               =       npsum;
ops.numpulsesperchunk       =       np;
ops.nChunks                 =       Nchunk;
ops.spkTh                   =       spkth;
plsout.ops                  =       ops;
plsout.chanMap              =       chanMap;
if ~isnan(stimpara.stimulus),        plsout.stimpara = stimpara; end

if ~exist(savingpath,'file')
    save(savingpath,'-v7.3','-struct','plsout');
end
%move(tsfolder,sn);

fprintf('Took %2.2f seconds...\n',toc);


end


function plotelectstim(chdatbefore, chdatafter, chtvec, meantempl, mt, plstempl ,ops, iter, savingpath)

savingpath = [savingpath,filesep, 'Chunk ', num2str(iter(2))];
if not(exist(savingpath, 'dir')); mkdir(savingpath); end
sfname = [num2str(iter(1),'%02g'),'-Artifact_removed_from_electrode_',num2str(iter(1))];
if exist([savingpath,filesep,sfname,'.png'],'file'), return; end

cols = lines(10);
tt = 0.03;
xax = linspace (0,ops.nt0 /ops.fs* 1e3, ops.nt0);
tvec = (double(chtvec) / ops.fs) / 600;

h = figure('pos',[10 200 1700 650],'color','w','vis','off');
subplot_tight(2,6,[1 5],tt);
plot(tvec, chdatbefore,'color',cols(6,:));
title(['Raw signal without stimulus subtraction for electrode: ',...
    num2str(iter(1)),', in chunck: ',num2str(iter(2))]);
% axis tight;
xlim([tvec(1), tvec(end)])
ax = gca;       ax.TickDir = 'out';     ax.TickLength = [0.005 0.005]; ax.XColor = 'none';
box off;        pbaspect([6 1 1]);


subplot_tight(2,6,[7 11],tt)
plot(tvec, chdatafter,'color',cols(7,:));
xticks(0:1:120);        xlabel('Time (min)');       box off;
set(gca,'YLim',ax.YLim, 'XLim',ax.XLim, 'TickDir','out','ticklength', [0.005 0.005]);
title(['Raw signal after  stimulus subtraction for electrode: ',...
    num2str(iter(1)),', in chunck: ',num2str(iter(2))]);
pbaspect([6 1 1]);


subplot_tight(2,6,6,tt);
plot(xax, plstempl','color',cols(6,:));         hold on;
plot(xax, mt','color','k','LineWidth',1);
set(gca,'YLim',ax.YLim, 'XLim',[0 max(xax)],'TickDir','out','ticklength', [0.025 0.025]);
box off;        axis square;        title('Stimulus pulses');


subplot_tight(2,6,12,tt);
plot(xax, transpose(plstempl-meantempl),'color',cols(7,:));
set(gca,'YLim',ax.YLim,'XLim',[0 max(xax)], 'TickDir','out','ticklength', [0.025 0.025]);
xticks(0:2:100);       box off;        axis square;     xlabel('Time (ms)');
title('Stimulus residuals');

%saveas(h,[savingpath,filesep, sfname,'.png']);
savepngFast(h, savingpath, sfname);
close all;
end




function [pulsetemplate, pulsestartstop, pulselocindex, besttemploc, spkth] = get_electrical_pulses(chdat, chtvec, ops)

[nChan, ~] = size(chdat);

redstep  = 5; % this is amount of reduction in threshold for each step
consitantpeaks = 0;     % this is to make sure that at least we find 5 consistant peaks
iter = 1;
npeaksfound = nan(1,50);

dataRAW =  (single(chdat')); % first in cpu then mover to gpu to avoid crashing
dataRAW = dataRAW./mad(dataRAW, 1, 1);
dataRAW = gpuArray(mean(dataRAW,2)); % take the average of all the signals across channels

while ~consitantpeaks
    
    smin = my_min(dataRAW, ops.loc_range, [1 2]);
    peaks = single(dataRAW<smin+1e-3 & dataRAW<ops.spkTh);
    
    if ops.spkTh > -15, redstep = 1; end
    ops.spkTh = ops.spkTh + redstep;
    if ops.spkTh > -2, break; end
    
    npeaksfound(iter) = double(gather(sum(peaks)));
    cp = npeaksfound(~isnan(npeaksfound));  cp = cp(cp>0);
    cpcounts = histc(cp,unique(cp)) > 2; %#ok % at least more that 2 time same number of peaks
    consitantpeaks = any(cpcounts);
    
    iter = iter+1;
    clearvars smin;
end
ops.spkTh = ops.spkTh - 5*redstep;
spkth = ops.spkTh;
p = gather(find(peaks)); % get all the peaks

if isempty(p) % in case there is no pulse, just return NaN
    [pulsetemplate, pulsestartstop, pulselocindex, besttemploc] = deal(nan(1,1)); 
    return;
end

pksdiff = diff(p);
p2toend = p(2:end);
pksbestloc = [ p(1,1); p2toend(pksdiff > ops.fs/10) ];
if isempty(pksbestloc), pksbestloc = NaN;  end
clearvars peaks dataRAW;

peakstimediff = (diff(round(diff(pksbestloc)/1000)*1000));
if not(all(peakstimediff==0))
    warning('inconsistant inter-peak-interval, some smelly shit is going on');
    disp(find(peakstimediff ~= 0));
end
besttemploc = gather(pksbestloc);

%np = [ops.minexpectedpulse-1, ops.minexpectedpulse+1]; % this is to dynamcally adapt the threshold
%redstep  = 10; % this is amount of reduction in threshold for each step
% while  any([length(unique(np)) > 1 , ~all(np > ops.minexpectedpulse)])
%     dataRAW = gpuArray (single(chdat'));
%     dataRAW = dataRAW./mad(dataRAW, 1, 1);
%     
%     smin = my_min(dataRAW, ops.loc_range, [1 2]);
%     peaks = single(dataRAW<smin+1e-3 & dataRAW<ops.spkTh);
%     
%     sigallranges = range(dataRAW,1);
%     sigstartranges = range(dataRAW(1:500,:),1);
%     sigtostimdiff =  gather(sigallranges - sigstartranges);
%     [~,goodchidx] = sort(sigtostimdiff,'descend');
%     
%     clearvars smin sigstartranges sigallranges dataRAW;
%     goodchpeaks = gather(peaks(:,goodchidx(1: ops.goodchnum))); % only check the top 5 channels
%     
%     pksbestloc = cell(1,ops.goodchnum);
%     for ii = 1:ops.goodchnum
%         p = find(goodchpeaks(:,ii));
%         if isempty(p), pksbestloc{ii} = NaN; continue; end  
%         pksdiff = diff(p);
%         p2toend = p(2:end);
%         pksbestloc{ii} = [ p(1,1); p2toend(pksdiff > ops.fs/10) ];
%     end
%     np = cellfun('length',pksbestloc);
%     if ops.spkTh > -15, redstep = 2; end
%     ops.spkTh = ops.spkTh + redstep;
%     if ops.spkTh > -6, break; end
% end
% spkth = ops.spkTh;
% besttemploc = median(cell2mat(pksbestloc),2);


correctedpeaks = false(size(chdat,2),size(chdat,1));
correctedpeaks( besttemploc, :) = true;
[row, col, ~] = find(correctedpeaks);

chdat = single (chdat);
clips = get_SpikeSample(chdat', row, col, ops, 0);

dt = 1:ops.nt0;
electwaveforms = squeeze(clips(:,:));

% temporal indices
indsT = repmat((row-ops.nt0min)', numel(dt), 1) + repmat(dt', 1, numel(row));
indsT = squeeze(indsT); % this is the indices of all the pulses
electpulsedat   = arrayfun(@(x)(electwaveforms(:,col==x)'),1:nChan,'un',0)';
% get time of each  electric event
% electpulsetimes = chtvec(indsT);
% pulsestartstop  = arrayfun(@(x)(electpulsetimes(:,col==x)'),1:nChan,'un',0)';

pulsestartstop = chtvec(indsT(:,col==1))'; % since it is same for all the channels
pulselocindex = indsT(:,col==1)';

pulsetemplate = cat(3, electpulsedat{:});

end
