

function InesHippocampusData(dp)

controldp = [dp,filesep,'01_control_chamber A\h5\'];
ectodp = [dp,filesep,'03_ectosomes_chamber C\h5\'];
exodp = [dp,filesep,'06_exosomes_chamber F\h5\'];

alldps = {controldp, ectodp, exodp};
alldata = cell(size(alldps));
for ii = 1:numel(alldps)
    alldata{ii} = load([alldps{ii},'ks_sorted\ksrasters\ksrasters.mat']);
    expnames = dir([alldps{ii},'*.h5']);
    alldata{ii}.expnames = {expnames.name}';
    alldata{ii}.datemodified = {expnames.date}';
    n = strrep(extractBetween( {expnames.name}','_','Hippo'),'T','-');
    alldata{ii}.expdate = datestr(datetime(n,'InputFormat','yyyy-MM-dd-HH-mm-ss'));
end

for ii = 1:numel(alldata)
    d = alldata{ii};
    dgoodidx = d.clusters(:,3)<5;
    fn = fieldnames(d);
    fn = fn(1:find(contains(fn,'template_info')));
    for jj = 1:numel(fn)
        if strcmp(fn{jj},'sort_params')
            d.sort_params.clusters_id = d.sort_params.clusters_id(dgoodidx,:);
        else
            d.(fn{jj}) = d.(fn{jj})(dgoodidx,:);
        end
    end
    alldata{ii} = d;
end



nbins = 100;
for ii = 1:numel(alldata)
    d = alldata{ii};
    d.nstim = size(d.spike_times,2);
    d.recduration = diff(d.sort_params.stim_start_end/d.sort_params.sampling_rate,1,2);
    d.meanamps = cellfun(@nanmean,d.amplitudes);
    d.frate = cellfun(@numel,d.spike_times) ./ d.recduration';
    d.numspikes = cellfun(@numel,d.spike_times);
    d.axonidx = contains([d.sort_info.comment],'axon');
    d.somaidx = ~d.axonidx;
    for jj = 1:d.nstim
        d.rasters{jj} =  CelltoMatUE(d.spike_times(:,jj))';
        bins = linspace(0,d.recduration(jj)+(d.recduration(jj)/nbins),nbins+1);
        p = nanmean(histc(d.rasters{jj},bins),2)'; %#ok
        d.psth(jj,:) = p(1:nbins);
        d.psthtimes(jj,:) = bins(1:nbins);
    end
    alldata{ii} = d;
end
%%
tt = 0.03;
for ii = 1:numel(alldata)
    d = alldata{ii};
    cols = lines(d.nstim);
    figure('pos',[200 0 1400 330*d.nstim],'color','w')
    for jj = 1:d.nstim
        subplot_tight(d.nstim*2,5,[(jj-1)*10+1, (jj-1)*10+4],[tt,tt])
        rasterPlotter(d.rasters{jj},[],cols(jj,:));
        set(gca,'xcolor','none','ticklength',[0.0025 0.0025]);
        box off;    %pbaspect([6 1 1]);
        ylabel('MEA electrodes');
        lb = unique([ceil(1:size(d.clusters,1)/4:size(d.clusters,1)),size(d.clusters,1)]);    yticks(lb);
        yticklabels(num2cellstr(d.clusters(lb,1)));
        axis([0 d.recduration(jj)+0.1 0 size(d.rasters{jj},2)+1]);
        title(strrep(d.expnames{jj},'_','-'));
        
        subplot_tight(d.nstim*2,5,[(jj-1)*10+6, (jj-1)*10+9],[tt+0.005, tt])
        plot(d.psthtimes(jj,:),d.psth(jj,:),'color',cols(jj,:));
        %pbaspect([6 1 1]);
        box off;    axis tight;
        set(gca,'ticklength',[0.0025 0.0025]);
        if jj == d.nstim, xlabel('time (sec)'); end
        ylabel('firing rate (Hz)');
        xlim([0 round(d.recduration(jj))]);
        xticks(round(0:d.recduration(jj)/2:d.recduration(jj)))
        ax = gca;   ax.Position(2)= ax.Position(2)+0.02;
        
    end
    subplot_tight(d.nstim*2,5,5,tt)
    thesis.colorboxplotOutline(d.frate,'color',cols);
    axis square;        box off;
    ylabel('firing rate (Hz)');
    
    subplot_tight(d.nstim*2,5,10,tt)
    cla
    thesis.colorboxplotOutline(d.frate,'color',cols);
    axis square;        box off;
    ylabel('amplitude (a.u.)');
    
    if d.nstim > 1
        subplot_tight(d.nstim*2,5,15,tt)
        thesis.colorboxplotOutline(d.numspikes,'color',cols);
        axis square;        box off;
        ylabel('spike numbers');
        xticks(1:d.nstim);
        
        subplot_tight(d.nstim*2,5,20,tt)
        line(1:d.nstim,d.meanamps,'marker','o','markersize',6,'linewidth',1,'color',[0 0 0 0.5],'markerfacecolor',0.5*[1 1 1])
        xlim([0.5 d.nstim+0.5])
        axis square;        box off;
        ylabel('amplitude (a.u.)');
        xticks(1:d.nstim);
    end
    
    
    if d.nstim > 2
        subplot_tight(d.nstim*2,5,25,tt)
        line(1:d.nstim,d.frate,'marker','o','markersize',6,'linewidth',1,'color',[0.8 0 0 0.5],'markerfacecolor',0.5*[1 0 0])
        xlim([0.5 d.nstim+0.5])
        axis square;        box off;
        ylabel('firing rate (Hz)');
        xticks(1:d.nstim);
        
        
        subplot_tight(d.nstim*2,5,30,tt)
        line(1:d.nstim,d.numspikes,'marker','o','markersize',6,'linewidth',1,'color',[0 0 0.8 0.5],'markerfacecolor',0.5*[0 0 1])
        xlim([0.5 d.nstim+0.5])
        axis square;        box off;
        ylabel('number of spikes');
        xticks(1:d.nstim);
        
    end
    
    
    saveas2(gcf,[dp,'\',cell2mat(extractBetween(alldps{ii},[dp,'\'],'\h5'))],600,'png');
    saveas2(gcf,[dp,'\',cell2mat(extractBetween(alldps{ii},[dp,'\'],'\h5'))],600,'pdf');
    close all;
    
end

%% population plots

lb = {'control','ectosome','exosome'};
figure('pos',[350 350 1500 500],'color','w')
subplot(1,3,1)
fr = {mean(alldata{1}.frate,2),mean(alldata{2}.frate,2),mean(alldata{3}.frate,2)};
thesis.colorboxplotOutline(fr);
axis square;    axis([0 4.5 -0.3 2]);
title('firing rate (Hz)');          xticklabels(lb);   xtickangle(45)
yticks(0:1:10);     box off;

subplot(1,3,2)
amp = {mean(alldata{1}.meanamps,2),mean(alldata{2}.meanamps,2),mean(alldata{3}.meanamps,2)};
thesis.colorboxplotOutline(amp);
axis square;    axis([0 4.5 -0.3 100]);
title('spike amplitude');          xticklabels(lb);   xtickangle(45)
yticks(0:50:1000);     box off;

subplot(1,3,3)
nspk = {mean(alldata{1}.numspikes,2),mean(alldata{2}.numspikes,2),mean(alldata{3}.numspikes,2)};
thesis.colorboxplotOutline(nspk);
axis square;    axis([0 4.5 -30 600]);
title('spike numbers');          xticklabels(lb);   xtickangle(45)
yticks(0:200:1000);     box off;

saveas2(gcf,[dp,'\population analysis'],600,'png');
saveas2(gcf,[dp,'\population analysis'],600,'pdf');
close all;

end