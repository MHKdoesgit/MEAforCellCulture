

function varargout = crossCorrAmpWaveformKSdata(datapath , kspath, stimid, varargin)
%
%%% crossCorrAmpWaveformKSdata %%%
%
%
% This function analyzes the spike trains of a selected stimulus and
% measure the corss correlation between all the channels and also auto
% correlations for each channel. Additionally, the function loads the spike
% waveforms of each unit from the kilosorted bin file and plots the
% electrical image, spike waveforms and spike amplitude in number of
% channels with strongest signal.
%
% ===============================Inputs====================================
%
%   datapath : path of the data folder.
%   kspath : path to ks_sorted forlder where alldata.bin file exists.
%   stimid : stimulus number to limit the analysis to one stimulus.
%
%================================Output====================================
%
%   acsqdata : array-structure containing auto/cross correlations, spike
%   waveforms, spikes amplitudes and spike templates.
%
% written by Mohammad, 06.05.2020.
%

totaltime = tic;
p = inputParser();      % check the user options.
p.addParameter('tbins', 0.5/1e3, @(x) isnumeric(x));
p.addParameter('lagnbins', 200, @(x) isnumeric(x));
p.addParameter('channelsaround', 15, @(x) isnumeric(x));
p.addParameter('elecImgAngle', 0, @(x) isnumeric(x));
p.addParameter('crosscorrthreshold', 0.25, @(x) isnumeric(x));
p.addParameter('plotscaling', 0.35, @(x) isnumeric(x));
p.addParameter('amplitudeNumPoints', 2500, @(x) isnumeric(x));
p.addParameter('waveformNumLines', 150, @(x) isnumeric(x));
p.parse(varargin{:});

rawdatpath = [datapath,filesep,'Data Analysis',filesep,'Raw Data',filesep];
stimnames = dir([rawdatpath,'*.mat']);      stimnames = {stimnames.name}';

% fist saving path
savingpath = [datapath,filesep,'Data Analysis',filesep,'Auto-Cross correlation and Sorting Quality Anaylsis',...
    filesep,extractBefore(stimnames{stimid},' for Exper'),filesep];
% folder making
if not(exist(savingpath,'dir')), mkdir(savingpath); end
if not(exist([savingpath,filesep,'acsq_data'],'dir')), mkdir([savingpath,filesep,'acsq_data']); end
% final file name
savefilename = [savingpath,filesep,'acsq_data',filesep,'autocorr spike quality for experiment on ',...
    datemaker(savingpath),'.mat'];
% if the analysis was done before, just plot the output
if exist(savefilename, 'file')
    analyzeflag = false;
    acsqdata = load(savefilename);
else
    analyzeflag = true;
end

if analyzeflag
    % updating parameters based on user input
    para.tbins = p.Results.tbins;
    para.lagnbins = p.Results.lagnbins;
    para.channelaround = p.Results.channelsaround;
    para.elecImgAngle = p.Results.crosscorrthreshold;
    para.ccthresh = p.Results.crosscorrthreshold;
    para.plotscaling = p.Results.plotscaling;
    para.ampNumPoints = p.Results.amplitudeNumPoints;
    para.wvformNumLines = p.Results.waveformNumLines;
    para.stimulusid = stimid;
    
    % loading ksrasters
    ksras = load([kspath,filesep,'ksrasters',filesep,'ksrasters.mat']);
    
    % load raw data
    rawdata = load([rawdatpath,stimnames{para.stimulusid}]);
    % loading clusters
    clus = rawdata.clusters.goodcells;
    
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
    % load raw data
    %rawdata = load([rawdatpath,stimnames{para.stimulusid}]);
    % auto and cross correlation
    ccdata = calcCCG(rawdata.spiketimes, para.tbins, para.lagnbins);
    % spikes waveforms
    spkwv = ksspkwaveforms(rawdata, kspath, ras.sort_params, ras.template_info, para.stimulusid, para.channelaround);
    % setting up the output
    acsqdata.ccdata = ccdata;
    acsqdata.spkwv = spkwv;
    acsqdata.ksrasters = ras;
    acsqdata.para = para;
    % saving
    save(savefilename,'-struct','acsqdata');
end
% plotting
plorACSQdata(acsqdata.ccdata, acsqdata.spkwv, acsqdata.ksrasters, acsqdata.para, savingpath);

disp(' And BOOM!!! Goes all the correlations');
disp(seconds2human (toc(totaltime)));
sound(struct2array(load('chirp.mat','y')));
% function output
varargout{1} = acsqdata;

end

%--------------------------------------------------------------------------------------------------%
%----------                               sub-functions                                  ----------%
%--------------------------------------------------------------------------------------------------%

function plorACSQdata(ccdata, spkwv, ras, para, savingpath)

tt= 0.05;
msg = [];
for ii = 1:size(ras.clusters,1)
    %
    h = figure('pos',[20 20 1800 900],'color','w','vis','off');
    % auto-correlation
    subplot_tight(3,3,1,tt)
    area(ccdata.laghalfms,ccdata.autocorrnopeak(ii,:),'FaceColor',rgb('dodgerblue'),'EdgeColor',rgb('royalblue'),'LineWidth',1)
    xlim([0 50]);       box  off;    xticks(0:25:100);
    title('auto-correlogram');      pbaspect([2 1 1]);      %xlabel('time (ms)');
    ax = gca;       ymax = ax.YLim(2);      ymax = ymax-ymax/10;
    chp = ras.sort_params.channel_positions;
    chp(:,2) = normalizetoRange(chp(:,2),ymax/2,ymax);
    chp(:,1) = normalizetoRange(chp(:,1),35,45);
    hold on;
    plot(chp(:,1),chp(:,2),'o','color','k','MarkerFaceColor',rgb('gray'),'MarkerSize',3);
    chid = ismember(ras.sort_params.channel_map,double(ras.clusters(ii,1)));
    plot(chp(chid,1),chp(chid,2),'o','color','r','MarkerFaceColor',rgb('red'),'MarkerSize',5);
    
    % spike template and electrical image
    subplot_tight(3,3,[2 5],tt)
    tmp = ras.template_info(ii);
    spktmpidx = spkwv.templateorder(ii,1:para.channelaround);
    coords = 1+squeeze(spkwv.channelcoords(ii,:,:))/1e2;
    if sum(spktmpidx) > 1
        tmp = tmp.templates(:,spkwv.templateorder(ii,1:para.channelaround));
        tmp = tmp(tmp(:,1)~=0,:);
        plotSpikeTemplate(squeeze(spkwv.templatesMean(ii,spkwv.channelsorder(ii,:),:))',coords,'k',para.elecImgAngle,'linewidth',1);
        hold on
        plotSpikeTemplate(tmp,coords,'r',para.elecImgAngle,'linewidth',1);
        arrdim = [ceil(sqrt(size(spkwv.coords,1))),ceil(sqrt(size(spkwv.coords,1)))];
        legend('mean response','template'); legend boxoff;      xticks(0:1:arrdim(1));     yticks(0:1:arrdim(2));
    end
    axis square;            ax = gca;
    
    % cross-correlation
    subplot_tight(3,3,4,tt)
    cc = squeeze(ccdata.crosscorr(ii,:,ccdata.lagamount+1:end));
    cc(ii,:) = NaN;
    ccviol = max(cc,[],2) > para.ccthresh;
    plot(ccdata.laghalfms,cc,'k');
    yAx = max(cc(:)+0.1);   if yAx <= 0 || isnan(yAx), yAx = 10; end
    axis([0 50 0 yAx]);
    title('cross-correlogram');
    if sum(ccviol)>1
        hold on;
        p = plot(ccdata.laghalfms,cc(ccviol,:),'r');
        legend(p,num2str(ras.clusters(ccviol,1:2))); legend boxoff;
        %title(['cross-correlogram (','violation: ',num2str(find(ccviol))',')']);
    end
    xlabel('time (ms)');        box off;        pbaspect([2 1 1]);      xticks(0:25:100);
    %clearvars cc ccviol p ax chp chid;
    
    % spike waveforms
    subplot_tight(3,3,[3 6],tt)
    swv = spkwv.spikeWaveforms{ii};
    xvals = coords(:,1); yvals = coords(:,2);
    alldist = sqrt( (xvals - xvals').^2 + (yvals - yvals').^2 );
    mindist = min(alldist(alldist(:)>0));
    if ~isempty(mindist) && size(swv,3) > 1
        tvec = linspace(-mindist*para.plotscaling, mindist*para.plotscaling, size(swv,2));
        ymax = max(abs(swv(1,:,:)),[],'all');
        xdata = xvals' + tvec';
        cols = gray(para.channelaround+4);
        if size(swv,3) > para.wvformNumLines, ns = para.wvformNumLines; else, ns = size(swv,3); end
        spkid = randperm(size(swv,3),ns);
        for jj = 1:para.channelaround
            ydat =  yvals(jj)' +  squeeze(swv(jj,:,spkid)) *(para.plotscaling)*mindist/ymax;
            x = repmat([xdata(:,jj);nan(1,1)],size(ydat,2),1);
            y= [ydat;nan(1,size(ydat,2))];
            plot(x(:),y(:),'color',cols(jj,:));
            %plot(xdata(:,jj),ydat,'color',cols(jj,:))
            hold on
        end
        %plotSpikeTemplate(squeeze(spkwv.templatesMean(ii,spkwv.channelsorder(ii,:),:))',coords,rgb('gold'),para.elecImgAngle,'linewidth',0.5);
        %clearvars swv xvals alldist mindist tvec ymax xdata spkid ns;
        axis square;         xticks(0:1:arrdim(1));     yticks(0:1:arrdim(2));      box off;
        xlim([ax.XLim]);     ylim([ax.YLim]);
    end
    
    % spike amplitudes
    subplot_tight(3,3,[7 9],tt)
    cols = repmat(0.5*[1 1 1],size(ras.amplitudes,2),1);%gray(size(ras.amplitudes,2));
    cols(para.stimulusid,:) = [1 0 0];
    st = ras.sort_params.stim_start_end;
    fs = ras.sort_params.sampling_rate*3600;
    for jj = 1:size(ras.amplitudes,2)
        amp = ras.amplitudes{ii,jj};
        
        if size(amp,1) > para.ampNumPoints, ns = para.ampNumPoints; else, ns = size(amp,1); end
        spkid = randperm(size(amp,1),ns);
        plot(linspace(st(jj,1),st(jj,2),ns)./fs,amp(spkid),...
            '.','MarkerEdgeColor',cols(jj,:),'MarkerFaceColor','none');
        hold on;
        xline(st(jj,1)/fs,'--r',num2str(jj,'%02d'));
        
    end
    ax = gca;       ax.YLim(1) = 2;             ax.TickLength = [0.001 0.001];
    ax.YTick = 0:4:24;      box off;
    xlabel('time (sec)');               ylabel('amplitude');
    
    % title and savenames
    [filename, filenamesavepng] = generateRGCname('Cross-Correlation, Auto Correlation and Sorting Quality', ras.sort_info(ii),savingpath);
    pngfilename = [num2str(ii,'%02d-'),extractAfter(filenamesavepng,'Cross-Correlation, Auto Correlation and ')];
    suptitle_mod(h,filename,2);
    % saving
    savepngFast(h,savingpath,pngfilename);
    close(h);
    fprintf(repmat('\b', 1, numel(msg)));
    msg = sprintf('Analysis for Cell %2d, Cluster %2d, is...wait for it...done!!!\n', ras.clusters(ii,1), ras.clusters(ii,2));
    fprintf(msg);
end

end

%--------------------------------------------------------------------------------------------------%

function res = ksspkwaveforms(rawdata, kspath, sortparams, kstemps, stimid, numchsig, varargin)

if nargin < 6, numchsig = 10;  end
tic;
Ncells = size(rawdata.clusters.goodcells,1);
binpath = fullfile(kspath, 'alldata.dat');

stimsamples = [sortparams.stim_start_end(2:end,1)-1;sortparams.stim_start_end(end,2)];
stimstarts = [0;stimsamples];
fs = sortparams.sampling_rate;

Tmin = ceil(-1e-3*fs); Tmax = floor(3*1e-3*fs);
dt = Tmin:Tmax;
Nt = numel(dt);
NchanTOT = sortparams.params_py.n_channels_dat;

Nspkmax = 2500 * 10;
Rmax = 10 * 60 * fs;
stimSums   = zeros(Ncells, NchanTOT, Nt, 1, 'single');
stimVars   = zeros(Ncells, NchanTOT, Nt, 1, 'single');
stimSpikes = zeros(Ncells,        1,  1, 1, 'single');
spikewaveforms = cell(Ncells,1);
channelsorder = zeros(Ncells, numchsig);
chlocations = zeros(Ncells, numchsig, 2);
spktemplateorder = zeros(Ncells, length(sortparams.channel_map));

%--------------------------------------------------------------------------
disp('Extracting electrical images...');
cmap = load(fullfile(kspath,'chanMap.mat'));
coords = [cmap.xcoords cmap.ycoords];
fid = fopen(binpath, 'r');

spiketimes = double(spikeCell2Mat(rawdata.spiketimes,fs));
spktimes = spiketimes(:,2);
spkids = spiketimes(:,1);
csamples = min(Rmax, stimsamples(stimid));

offset = 2 * NchanTOT*stimstarts(stimid);
fseek(fid, offset, 'bof');
dat = fread(fid, [NchanTOT csamples], '*int16');

spkids((spktimes+Tmax)>csamples | (spktimes+Tmin)<1) = [];
spktimes((spktimes+Tmax)>csamples | (spktimes+Tmin)<1) = [];

cspksall = accumarray(spkids, spktimes, [Ncells 1], @(x) {x});

for icell = 1:Ncells
    
    cellspikes = cspksall{icell};
    if isempty(cellspikes), continue; end
    [mx,m] = max(cellspikes);
    if mx+max(dt) > size(dat,2) % this is for rare cases where last spike is too close to the end of recording
        cellspikes = [cellspikes(1:m-1);cellspikes(m+1:end) ];
    end
    
    Nspikes = min(numel(cellspikes), Nspkmax);
    cellspikes = cellspikes(randperm(numel(cellspikes), Nspikes));
    
    spseek =  dt' + cellspikes';
    spkwvfrms = single(dat(:, spseek));
    spkwvfrms = reshape(spkwvfrms, NchanTOT, numel(dt), Nspikes);
    
    tmp = kstemps(icell);
    [~,tmporder] = sort(var(tmp.templates),'descend');
    [tf,chorder] = ismember(1:NchanTOT, sortparams.channel_map(tmporder(1:numchsig)));
    [~,p] = sort(chorder(tf));    idx = find(tf);    chorder = idx(p); % little trick to to unsort ismember output
    
    channelsorder(icell,:) = chorder;
    chlocations(icell, :, :) = coords(chorder,:);
    spktemplateorder(icell,:) = tmporder;
    
    spikewaveforms{icell} = spkwvfrms(chorder,:,:);
    stimSums  (icell, :, :, 1) = sum(spkwvfrms, 3);
    stimVars  (icell, :, :, 1) = var(spkwvfrms, [], 3);
    stimSpikes(icell, :, :, 1) = Nspikes;
end
fclose(fid);
%--------------------------------------------------------------------------
res.spikeWaveforms  = spikewaveforms;
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

