
function voltageTracesforPublication()

expdp = getExperimentPath('D:\Ines',1,false);

set(groot, 'defaultAxesTickDir', 'out');
set(groot,  'defaultAxesTickDirMode', 'manual');

stimnames = {'01_spontaneous_activity_5min';
    '02_cathodic-anodicpulse_200dur_100ppamp_50repeats';
    '03_cathodic-anodicpulse_200dur_160ppamp_50repeats';
    '04_cathodic-anodicpulse_200dur_200ppamp_50repeats';
    '05_cathodic-anodicpulse_200dur_400ppamp_50repeats';
    '06_cathodic-anodicpulse_200dur_600ppamp_50repeats';
    '07_cathodic-anodicpulse_200dur_800ppamp_50repeats';
    '08_cathodic-anodicpulse_200dur_1000ppamp_50repeats';
    '09_spontaneous_activity_5min'};
% 
% spondat = struct2array(load(expdp.matFilePaths{ismember(expdp.AnalysisFolderNames, 'Spontaneous_Activity_Analysis')}));
% esdat = struct2array(load(expdp.matFilePaths{ ismember(expdp.AnalysisFolderNames, 'Electric_Pulse_Burst_Analysis')}));
% 
% esdat = rmfield(esdat,{'interpulsebursts','pulsetime2peak'});

% ad = [spondat(1);esdat;spondat(end)];

selectidx = listdlg('PromptString','Select stimulus','Name','plot it, plot it baby!',...
    'SelectionMode','single','ListString',stimnames,'ListSize',[350,200],'OKString','Select');

% bd = ad(selectidx);
rawdat = dir(expdp.RawDataPaths{1});
rawdat = {rawdat(3:end).name}';
rawdat = load([expdp.RawDataPaths{1},filesep,rawdat{selectidx}]);
clus = rawdat.clusters.goodcells;

kspath = [expdp.exppaths{1},filesep,'h5/ks_sorted/'];
% loading ksrasters
ksras = load([kspath,filesep,'ksrasters',filesep,'ksrasters.mat']);
allids = [ksras.sort_info.id]';
goodids = ismember(allids,clus(:,4));
% clean up ksrasters
ras = ksras;
ras.amplitudes = ras.amplitudes(goodids,:);
ras.clusters = ras.clusters(goodids,:);
ras.sort_info = ras.sort_info(goodids,:);
ras.sort_params.clusters_id = ras.sort_params.clusters_id(goodids,:);
ras.spike_times = ras.spike_times(goodids,:);
ras.template_info = ras.template_info(goodids,:);
ras.sort_params = rmfield(ras.sort_params,{'similar_templates','templates_ind'});
ras.template_info = rmfield(ras.template_info,{'pc_feature_ind','template_feature_ind',...
    'cluster_group','cluster_quality','channel_id','template_ind_frq','template_order','channel_number','cluster_id'});
ras.charound  = 1;

spkwv = ksspkwaveforms(rawdat.spiketimes, kspath, ras.sort_params, ras.template_info, selectidx, ras.charound);


ellist = 11:88;
electrodelabels = compose('%d',reshape(ellist(not(mod(ellist,10)==0 | mod(ellist,10)==9 )),[8,8]));
analogchlabels = {'11','18','81','88'};
electrodelabels = electrodelabels(~ismember(electrodelabels,analogchlabels));

Ncells = size(clus,1);
Nchans = 60;
iter =  1;
lab = cell(Ncells,1);
for ii = 1:Nchans
    ch = spkwv.cellperchannels{ii};
    if ~isempty(ch)
        for jj = 1:numel(ch)
            lab{iter} = sprintf('%02d - Electrode %s, cell %g, cluster %d',iter, electrodelabels{ii},clus(iter,1),clus(iter,2));
            iter = iter +1;
        end
    end
end

cellidx = listdlg('PromptString','Select a cell for plotting','Name','Time for some data analysis!',...
    'SelectionMode','single','ListString',lab,'ListSize',[300,700],'OKString','Select');

p.tr = spkwv.traces(clus(cellidx,1),:);
p.tvec = spkwv.timevec;
p.idx = spkwv.cellperchannels{clus(cellidx,1)};
p.wav = spkwv.spikeWaveforms(p.idx);
p.wtvec = spkwv.spiketimevec(p.idx);

%%
xr = [0, max(p.tvec)];
yr = max(abs([floor(min(p.tr/10))*10, ceil(max(p.tr/10))*10 ]));
yr = [-yr, yr];
tt = [0.055, 0.04];
cols = lines(200); cols = cols(6:end,:);


h = figure('pos',[5 10 1850 950],'color','w','vis','on');

subplot_tight(3,1,1,tt)
plot(p.tvec, p.tr);
axis([xr,yr]);      %xticks(0:50:10000);     yticks(-1000:50:1000);
xlabel('Time (seconds)');       ylabel('Voltage (\muV)');
set(gca,'ticklength',[0.0025 0.0025]);        box off;
ax1 = gca;       ax1.Position(2) = ax1.Position(2)+0.05;

subplot_tight(3,1,2,tt)
plot(p.tvec, p.tr,'k');
hold on;
for jj = 1:numel(p.idx)
    plot(p.wtvec{jj}, p.wav{jj},'color',cols(jj,:));
end
axis([xr,yr]);      %xticks(0:50:10000);     yticks(-1000:50:1000);
xlabel('Time (seconds)');       ylabel('Voltage (\muV)');
set(gca,'ticklength',[0.0025 0.0025]);        box off;
ax2 = gca;       ax2.Position(2) = ax2.Position(2)+0.06;

for jj = 1:numel(p.idx)
    subplot_tight(3,numel(p.idx),numel(p.idx)*2+jj,tt)
    plot(linspace(0,4,size(p.wav{1},1)),p.wav{jj},'color',cols(jj,:));
    hold on
    plot(linspace(0,4,size(p.wav{1},1)),mean(p.wav{jj},2),'color','k','LineWidth',1);
    %axis([0,4,yr]);
    xlim([0 4]);
    xticks(0:2:100);     yticks(-1000:50:1000);
    title(lab{spkwv.cellperchannels{clus(cellidx,1)}(jj)});
    box off;
    pbaspect([2,1,1]);
end

filename = strrep([stimnames{selectidx},'Volrage traces for id ',extractBefore(lab{spkwv.cellperchannels{clus(cellidx,1)}(jj)},', cluster')...
    ' of ',expdp.expnumbers{1},' from experiment on ',expdp.expdates{1}],'_','-');
suptitle(h, filename,2);

isplotsaved = false;
savingpath = [expdp.mainfolder,filesep,'Population Analysis',filesep,'Voltage Traces',filesep];
if ~exist(savingpath, 'dir'), mkdir(savingpath); end
pngname = strrep(filename,',','_');     pngname = strrep(pngname,' ','_');

v = 0;
while ~isplotsaved
    anatype = questdlg('What?','save or what?','save','change axis','cancel','change axis');
    switch lower(anatype)
        case 'save'
            isplotsaved = true;
            savepngFast(h, savingpath, [pngname,num2str(v,'-v-%02d')]);
            saveas2([savingpath,filesep,[pngname,num2str(v,'-v-%02d')]],600,'pdf');
            close(h);
            close all;
            
        case 'change axis'
            definput = {[num2str(ax1.XLim(1)),' - ',num2str(ax1.XLim(2))],[num2str(ax1.YLim(1)),' - ',num2str(ax1.YLim(2))]};
            datrange = inputdlg({'Enter X-axis limits:','Enter Y-axis limits:'},'Select the range of your data',[1 35],definput);
            xr = strrep(datrange{1},' ','');   xr = sscanf(xr,'%f-%f');
            yr = strrep(datrange{2},' ','');   yr = sscanf(yr,'%f-%f');
            ax1.XLim = xr;      ax2.XLim = xr;
            ax1.YLim = yr;      ax2.YLim = yr;
            isplotsaved = false;
        case 'cancel'
            isplotsaved = true;
            close(h);
            close all; 
    end
    v = v+1;
end


end


%--------------------------------------------------------------------------------------------------%

function res = ksspkwaveforms(spktimesRaw, kspath, sortparams, kstemps, stimid, numchsig, varargin)

if nargin < 6, numchsig = 10;  end
tic;
Ncells = numel(spktimesRaw);
Nchans = sortparams.params_py.n_channels_dat;
binpath = fullfile(kspath, 'alldata.dat');

stimsamples = [sortparams.stim_start_end(2:end,1)-1;sortparams.stim_start_end(end,2)];
stimstarts = [0;stimsamples];
fs = sortparams.sampling_rate;

Tmin = ceil(-1e-3*fs); Tmax = floor(3*1e-3*fs);
dt = Tmin:Tmax;
Nt = numel(dt);
NchanTOT = sortparams.params_py.n_channels_dat;

Nspkmax = 250000 * 10;
Rmax = 20 * 60 * fs; % 20 min of data
stimSums   = zeros(Ncells, NchanTOT, Nt, 1, 'single');
stimVars   = zeros(Ncells, NchanTOT, Nt, 1, 'single');
stimSpikes = zeros(Ncells,        1,  1, 1, 'single');
[spikewaveforms, spktimevec] = deal(cell(Ncells,1));
channelsorder = zeros(Ncells, numchsig);
chlocations = zeros(Ncells, numchsig, 2);
spktemplateorder = zeros(Ncells, length(sortparams.channel_map));
chorderall = zeros(Nchans,1);
%--------------------------------------------------------------------------
disp('Extracting electrical images...');
cmap = load(fullfile(kspath,'chanMap.mat'));
coords = [cmap.xcoords cmap.ycoords];
fid = fopen(binpath, 'r');

spiketimes = double(spikeCell2Mat(spktimesRaw,fs));
spktimes = spiketimes(:,2);
spkids = spiketimes(:,1);
csamples = min(Rmax, (stimsamples(stimid)-stimstarts(stimid)));

offset = 2 * NchanTOT*stimstarts(stimid);
fseek(fid, offset, 'bof');
datraw = fread(fid, [NchanTOT csamples], '*int16');

% this is the conversion factor for the h5 msrd files.
% the value and code is taken from the  McsHDF5.McsAnalogStream lines 210
% to 251. These are all part of getConvertedData method.
conv_factor = 59605 /  1e6; % division to to convert voltage to microvolt
adzero = 0;
voltagedata = bsxfun(@minus,single(datraw),adzero); % this is only for MCS data manager
voltagedata = bsxfun(@times,voltagedata,conv_factor);

spkids((spktimes+Tmax)>csamples | (spktimes+Tmin)<1) = [];
spktimes((spktimes+Tmax)>csamples | (spktimes+Tmin)<1) = [];

cspksall = accumarray(spkids, spktimes, [Ncells 1], @(x) {x});

for icell = 1:Ncells
    
    cellspikes = cspksall{icell};
    if isempty(cellspikes), continue; end
    [mx,m] = max(cellspikes);
    if mx+max(dt) > size(voltagedata,2) % this is for rare cases where last spike is too close to the end of recording
        cellspikes = [cellspikes(1:m-1);cellspikes(m+1:end) ];
    end
    
    Nspikes = min(numel(cellspikes), Nspkmax);
    cellspikes = cellspikes(randperm(numel(cellspikes), Nspikes));
    
    spseek =  dt' + cellspikes';
    spkwvfrms = single(voltagedata(:, spseek));
    spkwvfrms = reshape(spkwvfrms, NchanTOT, numel(dt), Nspikes);
    
    tmp = kstemps(icell);
    [~,tmporder] = sort(var(tmp.templates),'descend');
    [tf,chorder] = ismember(1:NchanTOT, sortparams.channel_map(tmporder(1:numchsig)));
    [~,p] = sort(chorder(tf));    idx = find(tf);    chorder = idx(p); % little trick to to unsort ismember output
    
    channelsorder(icell,:) = chorder;
    chlocations(icell, :, :) = coords(chorder,:);
    spktemplateorder(icell,:) = tmporder;
    
    spikewaveforms{icell} = squeeze (spkwvfrms(chorder,:,:));
    spktimevec{icell} = spseek/fs;
    stimSums  (icell, :, :, 1) = sum(spkwvfrms, 3);
    stimVars  (icell, :, :, 1) = var(spkwvfrms, [], 3);
    stimSpikes(icell, :, :, 1) = Nspikes;
    chorderall(icell) = chorder;
end
fclose(fid);

cellsperch = cell(Nchans,1);
for ii = 1: Nchans
    cellsperch{ii} = find(chorderall == ii);
end

%--------------------------------------------------------------------------
tvec = (1:csamples)/fs;%linspace(stimstarts(stimid),stimsamples(stimid)-stimstarts(stimid),csamples);
res.traces = voltagedata;%(chorderall,:);
res.timevec = tvec;
res.spikeWaveforms  = spikewaveforms;
res.spiketimevec = spktimevec;
res.cellperchannels = cellsperch;
res.templatesMean   = stimSums./stimSpikes;
res.templatesStd    = sqrt(stimVars);
res.numSpikes      = stimSpikes;
res.templateTimes   = dt/fs;
res.channelsorder = channelsorder;
res.channelcoords = chlocations;
res.templateorder = spktemplateorder;
res.coords = coords*1e-6;
toc;
end
