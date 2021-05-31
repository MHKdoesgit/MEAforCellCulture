

function rastersPSTHforPublication()

ops.rasterXaxis = false;        % option to show x-axis for the raster plot
ops.rasterTitle = false;        % option to show the title for the raster plot
ops.BurstYellowShade = false;   % option to show a shaded region around detected pulses
ops.PSTHdetectedPeaks = false;  % option to show the detected peaks on the psth
ops.OnsetOffsetBurst = false;   % option to show the onset and offset of the detected bursts
ops.color = 'k';                % color of each plot. should be [r,g,b] and between 0 and 1

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

spondat = struct2array(load(expdp.matFilePaths{ismember(expdp.AnalysisFolderNames, 'Spontaneous_Activity_Analysis')}));
esdat = struct2array(load(expdp.matFilePaths{ ismember(expdp.AnalysisFolderNames, 'Electric_Pulse_Burst_Analysis')}));

esdat = rmfield(esdat,{'interpulsebursts','pulsetime2peak'});

ad = [spondat(1);esdat;spondat(end)];

selectidx = listdlg('PromptString','Select stimulus','Name','plot it, plot it baby!',...
    'SelectionMode','single','ListString',stimnames,'ListSize',[350,200],'OKString','Select');

bd = ad(selectidx);
rawdat = dir(expdp.RawDataPaths{1});
rawdat = {rawdat(3:end).name}';
rawdat = load([expdp.RawDataPaths{1},filesep,rawdat{selectidx}]);
clus = rawdat.clusters.goodcells;

anatype = questdlg('Select what to plot?','Rasters plot selection','All cells','Only bursty cells','All cells');

switch lower(anatype)
    case 'all cells'
        ras = bd.rastersall;
        psth = bd.psthall;
        psthtime = bd.psthtime;
        
    case 'only bursty cells'
        ras = bd.rastersbursty;
        psth = bd.psth;
        psthtime = bd.psthtime;
        clus = clus(bd.burstycellindex,:);
end

definput = {['0 - ',round(num2str(bd.psthtime(end)))],['0 - ',num2str(size(ras,1))]};
datrange = inputdlg({'Enter X-axis limits:','Enter Y-axis limits:'},'Select the range of your data',[1 35],definput);
xr = strrep(datrange{1},' ','');   xr = sscanf(xr,'%f-%f');
yr = strrep(datrange{2},' ','');   yr = sscanf(yr,'%f-%f');

% plotting goes here
if ops.BurstYellowShade
    pfun = @(bs,be,yax,varargin)(patch([bs,be,be,bs]',repmat([0 0 yax yax],size(bs,1),1)','y','edgecolor','none','facealpha',0.5));
end
tt = [0.055, 0.04];
%cols = lines(20); cols = cols(6:end,:);

h = figure('pos',[5 100 1875 850],'color','w','vis','on');
yAx = ceil(size(ras,1)/5)*5;
subplot_tight(2,1,1,tt)
if ops.BurstYellowShade
    pfun(bd.burststart,bd.burstend,yAx);        hold on;
end
rasterPlotter(ras',[],ops.color);             hold on;

if ops.OnsetOffsetBurst
    plot([bd.burststart bd.burststart],[0 yAx],'r--',[bd.burstend bd.burstend],[0 yAx],'g--');
end
%cellid = 12;
%rasterPlotter(ras{ii}(cellid,:)',[],'m',cellid);
if ops.rasterTitle
    titr = sprintf(['number of detected bursts: %d, peak threshold: %2.1f, burst per minute: %2.1f, mean firing rate: %2.1f (Hz),',...
        ' mean burst duration: %2.1f (ms), mean ibi: %2.1f, mean intra-burst frequency: %2.1f, reliablity Rsq: %2.2f'],...
        bd.numbursts,bd.peakthreshold, bd.indices.burstpermin,mean(bd.indices.firingRate),mean(bd.indices.burstduration),...
        mean(bd.indices.ibi), mean(bd.indices.intraburstfrq), bd.indices.Rsq);
    title(titr,'fontsize',9);
end
%axis([0 bd.psthtime(end) 0 yAx]);
axis([xr;yr(1);yr(2)+0.4]');
ylabel('Electrode nummbers');
if ops.rasterXaxis
    set(gca,'ticklength',[0.0025 0.0025]);
    xlabel('Time (seconds)');
    xticks(0:50:1000);
else
    set(gca,'ticklength',[0.0025 0.0025],'xcolor','none');
end
yt = ceil(0:size(ras,1)/4:size(ras,1));         yticks(yt);
yticklabels(clus([1,yt(2:end)],1)+1);   % pluse one matlab to python

subplot_tight(2,1,2,tt)
yAx = ceil(max(psth)/5)*5;
if ops.BurstYellowShade
    pfun(bd.burststart,bd.burstend,yAx);        hold on;
end
plot(psthtime, psth,'-','color',ops.color,'LineWidth',1);    %pbaspect([3,1,1]);
if ops.PSTHdetectedPeaks
    plot(bd.burstpeaktime,bd.burstpeaks,'ro','MarkerSize',4);
end

if ops.OnsetOffsetBurst
    plot([bd.burststart bd.burststart],[0 yAx],'r--',[bd.burstend bd.burstend],[0 yAx],'g--');
end
axis([0 bd.psthtime(end) 0 yAx]);   yticks(round(0:yAx/4:yAx));
xticks(0:50:1000);
xlabel('Time (seconds)');       ylabel('Firing rate (Hz)');
set(gca,'ticklength',[0.0025 0.0025]);        box off;
xlim(xr);
%title(titr,'fontsize',9);
filename = strrep([stimnames{selectidx},' rasters and PSTH for ',expdp.expnumbers{1},' from experiment on ',expdp.expdates{1}],'_','-');
suptitle(h, filename,2);

savingpath = [expdp.mainfolder,filesep,'Population Analysis',filesep,'Raster and PSTH plots',filesep];
if ~exist(savingpath, 'dir'), mkdir(savingpath); end

pngname = strrep(filename,',','_');     pngname = strrep(pngname,' ','_');
savepngFast(h, savingpath, pngname);
saveas2([savingpath,filesep,pngname],600,'pdf');
close(h);
close all;

end