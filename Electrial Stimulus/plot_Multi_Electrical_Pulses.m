

function plot_Multi_Electrical_Pulses(stimdat, ad, analysistype, savingpath, varargin)

set(groot, 'defaultAxesTickDir', 'out');
set(groot,  'defaultAxesTickDirMode', 'manual');

switch lower(analysistype)
    case 1
        plotburstanalysis(stimdat, ad, savingpath);
    case 2
        ploteachcellresponses(stimdat, ad, savingpath);
    case 3
        plotperpulseanalysis(stimdat, ad, savingpath)
end

end


function plotburstanalysis(stimdat, d, savingpath)

nstim = size(stimdat,1);
pfun = @(bs,be,yax,varargin)(patch([bs,be,be,bs]',repmat([0 0 yax yax],size(bs,1),1)','y','edgecolor','none','facealpha',0.5));
tt = 0.04;
cols = lines(20); cols = cols(6:end,:);
stimlist = cell(nstim,1);    for kk= 1:nstim, stimlist{kk} = stimdat(kk).stimPara.stimulus; end
%msg = [];

if ~exist(savingpath,'dir'), mkdir(savingpath); end

for ampid = 1: nstim
    
    %clus = stimdat(ampid).clusters.goodcells;
    %para = stimdat(ampid).stimPara;
    %chinfo = stimdat(ampid).sortinginfo;
    expdate = stimdat(ampid).date;
    stim = strrep(stimdat(ampid).stimulus,'_',',');     stim(3)='-';
    
    bd = d(ampid);
    
    ras = {bd.rastersall, bd.rastersbursty};
    pboth = {bd.psthall, bd.psth};
    expfn = {'all cells','only bursty cells'};
    
    
    for ii = 1:numel(ras)
        
        h = figure('pos',[5 100 1875 850],'color','w','vis','off');
        yAx = size(ras{ii},1)+3;
        subplot_tight(2,1,1,tt)
        pfun(bd.burststart,bd.burstend,yAx);        hold on;
        rasterPlotter(ras{ii}',[],'k');             hold on;
        %plot([bd.burststart bd.burststart],[0 yAx],'r--',[bd.burstend bd.burstend],[0 yAx],'g--');
        %cellid = 12;
        %rasterPlotter(ras{ii}(cellid,:)',[],'m',cellid);
        titr = sprintf(['number of detected bursts: %d, peak threshold: %2.1f, burst per minute: %2.1f, mean firing rate: %2.1f (Hz),',...
            ' mean burst duration: %2.1f (ms), mean ibi: %2.1f, mean intra-burst frequency: %2.1f, reliablity Rsq: %2.2f'],...
            bd.numbursts,bd.peakthreshold, bd.indices.burstpermin,mean(bd.indices.firingRate),mean(bd.indices.burstduration),...
            mean(bd.indices.ibi), mean(bd.indices.intraburstfrq), bd.indices.Rsq);
        title(titr,'fontsize',9);
        axis([0 bd.psthtime(end) 0 yAx]);   ylabel('number of cells');
        set(gca,'ticklength',[0.0025 0.0025],'xcolor','none');
        
        
        subplot_tight(2,1,2,tt)
        yAx = ceil(max(pboth{ii})/5)*5;
        pfun(bd.burststart,bd.burstend,yAx);        hold on;
        plot(bd.psthtime, pboth{ii},'k-');    %pbaspect([3,1,1]);
        if ii ==2
            plot(bd.burstpeaktime,bd.burstpeaks,'ro','MarkerSize',4);
        end
        %plot([bd.burststart bd.burststart],[0 yAx],'r--',[bd.burstend bd.burstend],[0 yAx],'g--');
        axis([0 bd.psthtime(end) 0 yAx]);   yticks(round(0:yAx/4:yAx));
        xlabel('time (seconds)');       ylabel('Firing rate (Hz)');
        set(gca,'ticklength',[0.0025 0.0025]);        box off;
        %title(titr,'fontsize',9);
        filename = [stim,' burst analysis for ',expfn{ii},' from experiment on ',expdate];
        suptitle(h, filename,2);
        
        pngname = strrep(filename,',','_');     pngname = strrep(pngname,' ','_');
        savepngFast(h, savingpath, pngname);
        close(h);
        
    end
    
end

%---------------------------------------------------------------------------------------------------


idx = [d.indices]';
fn = fieldnames(idx);
fn = fn(~ismember(fn,{'pearsoncoeff','explainedvar','Rsq_allcells','pearsoncoeff_allcells','explainedvar_allcells'}));
for jj = 1:numel(fn)
    popidx.(fn{jj}) =  CelltoMatUE({idx.(fn{jj})});
end
popidx.numbursts = [d.numbursts];
popidx.peakthreshold = [d.peakthreshold];
%popidx.spkamps = CelltoMatUE({d.meanamplitudes});
fn = fieldnames(popidx);

if all(strcmpi(stimlist,'spontaneous'))
    stimamps = 1:nstim;
    xaxticklab = {'before','after'};
    xaxlab = [];
    titr = ['Population analysis for spontaneous acitivty of all the recorded cells from in experiment at ',expdate];
else
    stimamps = [stimdat.stimPara];
    stimamps = [stimamps.amplitude];
    xaxticklab = cellstr(num2str(transpose(stimamps)))';
    xaxlab = 'Stimulus amplitudes';
    titr = ['Population analysis for electrical stimulation of all the recorded cells from in experiment at ',expdate];
    
end

h = figure('pos',[10 50 1500 950],'color','w','vis','off');
for jj = 1:16
    subplot_tight(4,4,jj,0.05)
    p = popidx.(fn{jj});
    switch lower(fn{jj})
        case 'burstpermin', lgtxt = 'Bursts/minute';
        case 'firingrate',  lgtxt = 'Firing rate (Hz)';
        case 'ibi',  lgtxt = 'Inter-burst-intervals (sec)';
        case 'intraburstfrq',  lgtxt = 'Intera-burst spiking frequency (Hz)';
        case 'spkinburstpercent',  lgtxt = 'Spikes in burst (%)';
        case 'burstduration',  lgtxt = 'Burst duration (ms)';
        case 'numbursts',  lgtxt = 'Total detected bursts';
        case 'peakthreshold', lgtxt = 'PSTH dectection threshold';
        case 'rsq',  lgtxt = 'Reliability R^2';
        case 'totalspikes',   lgtxt = 'Total number of spikes (all cells)';
        case 'totalburstyspikes',   lgtxt = 'Total number of spikes (bursty cells)';
        case 'allchanfiringrate',   lgtxt = 'Firing rate (all cells)';
        case 'allmeanisi',   lgtxt = 'Average inter-spike-intervals';
        case 'meanisibursty',   lgtxt = 'Average inter-spike-intervals (bursty cells)';
        case 'allchansmeanamplitudes',   lgtxt = 'Average spike amplitudes (all cells)';
        case 'meanamplitudes',   lgtxt = 'Average spike amplitudes (bursty  cells)';
        otherwise, lgtxt = fn{jj};
    end
    
    if size(p,1) == 1
        ci = randi(10);
        plot(stimamps,p,'-o','MarkerFaceColor',cols(ci,:),'Color',cols(ci,:),'LineWidth',2);
        ax = gca;        xticks(stimamps);
        xlim([min(stimamps)-mean(diff(stimamps))/4, max(stimamps)+mean(diff(stimamps))/4]);
        if strcmpi(fn{jj},'burstpermin'), ax.YLim(1) = 0; end
    else
        datapointboxplot(p,'Colors',cols,'BoxWidth',0.3);
        ax = gca;
        if strcmpi(fn{jj},'intraburstfrq') && max(ax.YLim) > 300, ylim([-5 300]); end
        if strcmpi(fn{jj},'burstduration') && max(ax.YLim) > 800, ylim([-5 300]); end
        %if strcmpi(fn{jj},'ibi') && max(ax.YLim) > 800, ylim([-5 300]); end
    end
    lgspaces = strfind(lgtxt,' ');
    if any(lgspaces > 24), lgtxt = {lgtxt(1:(lgspaces(ceil(length(lgspaces)/2))-1)),lgtxt(lgspaces(ceil(length(lgspaces)/2))+1:end)}; end
    xlabel(xaxlab);       xticklabels(xaxticklab);
    box off;             ylabel(lgtxt);     pbaspect([1.2,1,1]);
    ax.XTickLabelRotation = 45;
end
suptitle(h, titr,2);

savepngFast(h, savingpath, titr);
close(h);

%---------------------------------------------------------------------------------------------------

if all(strcmpi(stimlist,'spontaneous'))
    nrows = 1;
    ncols = 3;
    titr = ['PSTHs of all bursts from spontaneous acitivty of all the recorded cells in experiment at ',expdate];
else
    nrows = 2;
    ncols = 4;
    titr = ['PSTHs of all bursts from electrical stimulation of all the recorded cells in experiment at ',expdate];
end

h = figure('pos',[50 100 1750 nrows*500],'color','w','vis','off');
yaxall = zeros(1,nstim);
for ampid = 1: nstim
    
    stim = strrep(stimdat(ampid).stimulus,'_',',');     stim(3)='-';
    
    bd = d(ampid).raspsthallbursts;
    subplot_tight(nrows,ncols,ampid,tt)
    yax = ceil(max(bd.psth)/10)*10;
    xax = ceil(max(bd.psthtime));
    bar(bd.psthtime,bd.psth,'BarWidth',1,'FaceColor',cols(ampid,:),'EdgeColor',abs(cols(ampid,:)-0.2));
    hold on;
    plot(bd.decaytvec/1e3,bd.decayval,'k','LineWidth',1);
    plot(bd.latency/1e3,bd.peak,'ro','LineWidth',2);
    text(xax-0.5,yax-(yax/4),{['Decay index: ',...
        num2str(round(bd.psthdecayindex,1))],['Peak: ',num2str(round(bd.peak,1))]})
    axis([0 xax 0 yax]);        xlabel('Time(sec)');        ylabel('Firing rate (Hz)');
    box off;
    xticks(0:0.2:10);        yticks(ceil(0:yax/4:yax));     pbaspect([2,1,1]);
    title({stim,[' ,Rsq: ',num2str(round(bd.rsq,2)),', pearson-coeff: ',num2str(round(bd.pearsoncoeff,2))]},'FontSize',9)
    yaxall(ampid) = yax;
end
subplot_tight(nrows,ncols,ampid+1,tt)
lgtxt = cell(1,nstim);
for ampid = 1: nstim
    plot(d(ampid).raspsthallbursts.psthtime,d(ampid).raspsthallbursts.psth,'Color',cols(ampid,:),'LineWidth',1.5);     hold on;
    lgtxt{ampid} = ['Stimulus ',num2str(ampid,'%02g')];
end
yax = max(yaxall);
axis([0 xax 0 yax]);        xlabel('Time(sec)');        ylabel('Firing rate (Hz)');
xticks(0:0.2:10);        yticks(ceil(0:yax/4:yax));     pbaspect([2,1,1]);
legend(lgtxt);      legend boxoff;          box off;
title('All psths together','FontSize',9)

suptitle(h, titr,2);
savepngFast(h, savingpath, titr);
close(h);




end

function ploteachcellresponses(stimdat, d, savingpath)

ncells = size(d.rasters,1);
nstim = size(stimdat,1);
clus = stimdat(1).clusters.goodcells;
sortinfo = stimdat(1).sortinginfo;
expdate = stimdat(1).date;

ncols = 3;
nrows = 7;
tt = [0.025 0.05];
cols = lines(12); cols = cols(6:end,:);
stimamps = [stimdat.stimPara];      stimamps = [stimamps.amplitude];
msg = [];


for ii = 1: ncells
    
    h = figure('pos',[50 50 1700 900], 'color','w','vis','off');
    
    yax = ceil(max(max(d.psth(ii,:,:),[],2))/2)*2;
    if yax <1, yax = 1; end
    for ampid =1:  nstim
        
        electstim = stimdat(ampid).electricalstimulus;
        numstim = numel(electstim.pulsedetectedat);
        stimdur = d.stimdur;
        tb = d.timebeforespk;
        stim = strrep(stimdat(ampid).stimulus,'_',',');     stim(3)='-';
        
        subplot_tight(nrows, ncols, ncols*(ampid-1)+1, tt )
        rasterPlotter(d.beforerasters{ii,ampid}',[],0.4*[1 1 1])
        rasterPlotter(d.timebeforespk + d.rasters{ii,ampid}',[],cols(2,:))
        xline(d.timebeforespk,'-.b');
        axis([0 stimdur+tb 0 numstim+1]);
        yticks(0:25:100);        ylabel('Trials');
        %if ampid==1, title(['Cell ID: ', num2str(ii)]); end
        if ampid == nstim
            xticks(0:1:numstim);     xticklabels(-tb:1:numstim-tb);
            xlabel('Time(sec)');
        else
            set(gca,'xcolor','none');
        end
        set(gca,'tickdir','out');
        
        subplot_tight(nrows, ncols, ncols*(ampid-1)+2, tt )
        plot(d.beforepsthtime,squeeze(d.beforepsth(ii,:,ampid)),'color',0.4*[1 1 1]);
        hold on;
        plot(d.psthtime + tb,squeeze(d.psth(ii,:,ampid)),'color',cols(2,:),'linewidth',1);
        axis([0 stimdur+tb 0 yax]);
        xline(d.timebeforespk,'-.b');
        box off;
        set(gca,'tickdir','out');
        if ampid == nstim
            xticks(0:1:numstim);     xticklabels(-d.timebeforespk:1:numstim-d.timebeforespk);
            xlabel('Time(sec)');
        else
            set(gca,'xcolor','none');
        end
        title(stim,'fontsize',8);
    end
    subplot_tight(3, ncols, 3, 2*tt )
    plot(stimamps,d.peak(ii,:),'-o','MarkerFaceColor',cols(3,:),'Color',abs(cols(3,:)-0.1),'LineWidth',1.5);
    box off;        title('Peak');     ylabel('Firing rate (Hz)');
    
    subplot_tight(3, ncols, 6, 2*tt )
    plot(stimamps,d.latency(ii,:),'-o','MarkerFaceColor',cols(4,:),'Color',abs(cols(4,:)-0.1),'LineWidth',1.5);
    box off;        title('Time-to-peak (ms)');     ylabel('Time (ms)');
    subplot_tight(3, ncols, 9, 2*tt )
    plot(stimamps,d.decaytime(ii,:),'-o','MarkerFaceColor',cols(5,:),'Color',abs(cols(5,:)-0.1),'LineWidth',1.5);
    box off;        title('Decay time');     ylabel('Decay-coeff');
    ax = gca;   if max(ax.YLim) > 500, ax.YLim = [0 500]; end
    xlabel('Stimulus amplitude');
    
    [filename,pngname] = cell_label('Electrical Stimulus', sortinfo(ii), expdate, ii);
    suptitle(h, filename,2);
    savepngFast(h, savingpath, pngname);
    %----------------------------------------------------------------------
    fprintf(repmat('\b', 1, numel(msg)));
    msg = sprintf('Cell %d/%d. Time elapsed %2.2f s...\n', ii,size(clus,1),toc);
    fprintf(msg);
    close(h);
    
end


end

function plotperpulseanalysis(stimdat, d, savingpath)

nstim = size(stimdat,1);
expdate = stimdat(1).date;
tt = 0.01;
cols = lines(24); cols = cols(6:end,:);
msg = [];

for ampid = 1:nstim
    
    clus = stimdat(ampid).clusters.goodcells;
    para = stimdat(ampid).stimPara;
    %chinfo = stimdat(ampid).sortinginfo;
    
    stim = strrep(stimdat(ampid).stimulus,'_',',');     stim(3)='-';
    filename = [stim,' rasters showing pulse effects for experiment on ', expdate];
    pngname = strrep(filename,',','_');     pngname = strrep(pngname,' ','_');
    
    ncols = 10;
    nrows = para.repeats/ncols;
    %---------------------------------------------------------------------------------------------------
    h = figure('pos',[1 30 1900 970], 'color','w','vis','off');
    for ii = 1:para.repeats
        subplot_tight(nrows,ncols,ii,tt)
        rasterPlotter(d(ampid).rasters{ii},[],cols(2,:));
        %axis square ;
        pbaspect([1.1 1 1]);
        ylim([0 ceil(size(clus,1)/2)*2]);
        title(ii);
        if mod(ii,ncols)==1, set(gca,'xtick',[],'ytick',0:size(clus,1)/2:size(clus,1));
        else,            set(gca,'xtick',[],'ytick',[]);
        end
        if ii > para.repeats - ncols, set(gca,'xtick',0:2:20); end
        box on;
    end
    suptitle(h, filename,2);
    savepngFast(h, savingpath, pngname);
    close(h);
    %---------------------------------------------------------------------------------------------------
    h = figure('pos',[1 30 1900 970], 'color','w','vis','off');
    yax = ceil(max(d(ampid).psth,[],'all')/2)*2;
    for ii = 1: para.repeats
        subplot_tight(5,10,ii,0.01)
        plot(d(ampid).psthtime,d(ampid).psth(ii,:),'color','k','LineWidth',1);
        pbaspect([1.1 1 1]);    title(ii);
        axis([0 round(d(ampid).stimdur) 0 yax+2]);
        if mod(ii,10)==1,        set(gca,'ytick',0:yax/2:yax);   else, set(gca,'ytick',[]);    end
        if ii > 40, set(gca,'xtick',0:5:100); else, set(gca,'xtick',[]); end
        box on
        
    end
    filename = [stim,' psths showing pulse effects for experiment on ', expdate];
    suptitle(h, filename,2);
    pngname = strrep(filename,',','_');     pngname = strrep(pngname,' ','_');
    savepngFast(h, savingpath, pngname);
    %----------------------------------------------------------------------
    fprintf(repmat('\b', 1, numel(msg)));
    msg = sprintf('Analysis for experiment %d is finito.... Time elapsed %2.2f s...\n', ampid,toc);
    fprintf(msg);
    close(h);
    
    
end

%---------------------------------------------------------------------------------------------------
ras = [{d.allrasters};{d.burstyrastes}];
xax = min([d.recordingdur]);
[yax,~] = cellfun(@size,ras);
yax = ceil(yax/2)*2;
titr = {'all recorded cells','only bursty cells'};

for ii = 1:2
    h = figure('pos',[50 50 1700 900],'color','w','vis','off');
    for ampid = 1: nstim
        subplot_tight(nstim,1,ampid,0.025)
        stim = strrep(stimdat(ampid).stimulus,'_',',');     stim(3)='-';
        rasterPlotter(ras{ii,ampid}',[],'k');       hold on;
        plot([d(ampid).electstimtime, d(ampid).electstimtime],[0 yax(ii,ampid)+2],'r-');
        axis([0 xax 0 yax(ii,ampid)+2]);
        box off;        set(gca,'xcolor','none','ytick',0:yax(ii,ampid)/2:yax(ii,ampid),'ticklength',[0.002 0.002]);
        if ampid == nstim, set(gca,'xcolor','k'); end
        ylabel('# cells');
        if ii==2
            title([stim,' (', num2str(yax(2,ampid)),' out of ',num2str(yax(1,ampid)),' cells)'],'fontsize',8);
        else
            title(stim,'fontsize',8);
        end
    end
    filename = ['Rasters for ',titr{ii},' from expreiment on ',expdate];
    suptitle(h,filename,2);
    pngname = strrep(filename,',','_');     pngname = strrep(pngname,' ','_');
    savepngFast(h, savingpath, pngname);
    close(h);
end


%---------------------------------------------------------------------------------------------------
p = [{d.allpsth}; {d.burstypsth}];

xax = min([d.recordingdur]);
[yax,~] = cellfun(@size,ras);
[yaxm,~] = cellfun(@max,p);
yaxm = ceil(yaxm/2)*2;
titr = {'all recorded cells','only bursty cells'};

for ii = 1:2
    h = figure('pos',[50 50 1700 900],'color','w','vis','off');
    for ampid = 1: nstim
        subplot_tight(nstim,1,ampid,0.025)
        stim = strrep(stimdat(ampid).stimulus,'_',',');     stim(3)='-';
        tvec = d(ampid).allpsthtime;
        plot(tvec,p{ii,ampid},'k');     hold on;
        %rasterPlotter(ras{ii,ampid}',[],'k');       hold on;
        plot([d(ampid).electstimtime, d(ampid).electstimtime],[0 yaxm(ii,ampid)+2],'r-');
        axis([0 xax 0 yaxm(ii,ampid)+2]);
        box off;        set(gca,'xcolor','none','ytick',0:yaxm(ii,ampid)/2:yaxm(ii,ampid),'ticklength',[0.002 0.002]);
        if ampid == nstim, set(gca,'xcolor','k'); end
        ylabel('Rate (Hz)');
        if ii==2
            title([stim,' (', num2str(yax(2,ampid)),' out of ',num2str(yax(1,ampid)),' cells)'],'fontsize',8);
        else
            title(stim,'fontsize',8);
        end
    end
    filename = ['PSTHs for ',titr{ii},' from expreiment on ',expdate];
    suptitle(h,filename,2);
    pngname = strrep(filename,',','_');     pngname = strrep(pngname,' ','_');
    savepngFast(h, savingpath, pngname);
    close(h);
end

%---------------------------------------------------------------------------------------------------

for ampid = 1: nstim
    es = d(ampid).esstim;
    npulses = size(es.pulsedetectedat,1);
    yax = max(es.pulsetemplates,[],'all');
    xax = size(es.pulsetemplates,2);
    stim = strrep(stimdat(ampid).stimulus,'_',',');     stim(3)='-';
    
    h = figure('pos',[150 50 1500 900],'color','w','vis','off');
    
    for ii = 1: npulses
        subplot_tight(5,10,ii,0.01)
        plot(squeeze(es.pulsetemplates(ii,:,:)))
        axis([0 xax -yax yax]);
        axis square off;
        title(['pulse ',num2str(ii)])
    end
    filename = [stim,' electrical pulses for expreiment on ',expdate];
    suptitle(h,filename,2);
    pngname = strrep(filename,',','_');     pngname = strrep(pngname,' ','_');
    savepngFast(h, savingpath, pngname);
    close(h);
    
    %---------------------------------------------------------------------------------------------------
    essub = es.pulsetemplates - es.mean2subtract;
    yax = max(abs(essub),[],'all');
    xax = size(essub,2);
    h = figure('pos',[150 50 1500 900],'color','w','vis','off');
    
    for ii = 1: npulses
        subplot_tight(5,10,ii,0.01)
        plot(squeeze(essub(ii,:,:)))
        axis([0 xax -yax yax]);
        axis square off;
        title(['pulse ',num2str(ii)])
    end
    filename = [stim,' residuals for expreiment on ',expdate];
    suptitle(h,filename,2);
    pngname = strrep(filename,',','_');     pngname = strrep(pngname,' ','_');
    savepngFast(h, savingpath, pngname);
    close(h);
end

%---------------------------------------------------------------------------------------------------
h = figure('pos',[50 50 1700 900],'color','w','vis','off');
for ampid = 1: nstim
    subplot_tight(nstim,1,ampid,0.03)
    stim = strrep(stimdat(ampid).stimulus,'_',',');     stim(3)='-';
    nump = length(d(ampid).pulsetime2peak);
    stem(1:nump, d(ampid).pulsetime2peak,'color',cols(ampid,:),'markerfacecolor',cols(ampid,:),'LineWidth',1)
    axis([0 nump 0 5]);
    box off;        set(gca,'xcolor','none','ytick',0:1:10,'ticklength',[0.002 0.002]);
    if ampid == nstim, set(gca,'xcolor','k'); end
    ylabel('Delay (sec)');
    title(stim);
    box off;
end
filename = ['Time delay between each pulse and the following burst from expreiment on ',expdate];
suptitle(h,filename,2);
pngname = strrep(filename,',','_');     pngname = strrep(pngname,' ','_');
savepngFast(h, savingpath, pngname);
close(h);


end