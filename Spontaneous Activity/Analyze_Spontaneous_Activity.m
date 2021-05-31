

function varargout = Analyze_Spontaneous_Activity(datapath, varargin)



totaltime =tic;
if (nargin < 1),    datapath = uigetdir();  end


rawdatpath = [datapath, filesep, 'Data Analysis/Raw Data', filesep];
stimlist = dir( [rawdatpath, '*spont*.mat']);
stimlist = {stimlist.name}';
stimnums = str2double(extractBefore(stimlist,'_spon'));

stimdat = cell(size(stimlist));
for ii = 1: numel(stimlist)
    stimdat{ii} = load([rawdatpath, filesep, stimlist{ii}]);
end
stimdat = cell2mat(stimdat);    

recordingdur = diff(stimdat(1).stim_start_end,1,2) ./ stimdat(1).samplingrate;
recordingdur = recordingdur(stimnums);
[burstactivity]= deal(cell(numel(stimlist),1));

for ii = 1: numel(stimlist)
    thisExp = stimdat(ii);
    burstactivity{ii} = burst_analysis(thisExp.spiketimes,recordingdur(ii));
    amps = stimdat(ii).amplitudes;
    burstactivity{ii}.amplitudes = amps(burstactivity{ii}.burstycellindex);
    burstactivity{ii}.allChansAmplitudes = amps;
    burstactivity{ii}.indices.meanAmplitudes = cellfun(@mean,amps(burstactivity{ii}.burstycellindex));
    burstactivity{ii}.indices.allChansMeanAmplitudes = cellfun(@mean,amps);
    burstactivity{ii}.spiketimes = stimdat(ii).spiketimes(burstactivity{ii}.burstycellindex);
end
burstactivity = cell2mat(burstactivity);

savingpath = [datapath, filesep, 'Data Analysis',filesep,'Spontaneous_Activity_Analysis'];
if ~exist(savingpath,'dir'), mkdir(savingpath); end

filename = ['Spontaneous activity analysis for experiment on ', stimdat(1).date];
save([savingpath,'\',filename,'.mat'],'-v7.3','burstactivity');

plot_Multi_Electrical_Pulses(stimdat, burstactivity,1,savingpath);

sound(struct2array(load('gong.mat','y')));
disp(seconds2human (toc(totaltime)));
varargout{1} = burstactivity;

% [thisExp, savingpath] = loadRawData(datapath,'spon','Spontaneous_Activity_Analysis','spondata');
% 
% allrasters = CelltoMatUE(thisExp.spiketimes);
% binlength = 100/1e3;
% 
% clus = thisExp.clusters.goodcells;
% endofexp = ceil(nanmax(allrasters,[],'all')/10)*10;
% 
% binVec = linspace(0,endofexp,1+endofexp/binlength);
% psth = mean(histc(allrasters',binVec),2)*(1/binlength); %#ok
% 
% chidx = find (diff(clus(:,1))>0)+1;
% %%
% h = figure('pos',[200 200 1920 1080],'color','w','vis','on');
% subplot_tight(2,1,1,0.05)
% rasterPlotter(allrasters',[],'k');
% if size(clus,1) > 100
%     hold on;
%     plot([0 endofexp],[chidx chidx],'r--');
% end
% ax = gca;
% ax.XColor = 'none';
% ax.TickDir = 'out';     ax.TickLength = [0.005 0.05];
% axis([0 endofexp 0 size(clus,1)+1]);
% %yticks(0 : size(clus,1)/6 : size(clus,1))
% yticks(chidx(1:floor(length(chidx)/6):end));
% yticklabels (clus(chidx(1:floor(length(chidx)/6):end),1));
% ylabel('number of electrodes');
% title('rasters');
% 
% subplot_tight(2,1,2,0.05)
% plot(binVec/60,psth,'color',rgb('dodgerblue'),'linewidth',1);
% ax = gca;       box off;
% ax.TickDir = 'out';     ax.TickLength = [0.005 0.05];
% yAx = ceil(max(psth)/2)*2;
% axis([0 endofexp/60 0 yAx]);
% 
% xticks(0:1:100);
% xlabel('time (min)');
% yticks(0:yAx/2:yAx);
% ylabel('rate (Hz)');
% title(['psth (',num2str(binlength*1e3),' ms bins)']);
% 
% %[filename, filenamesavepng] = generateRGCname('Cross-Correlation, Auto Correlation and Sorting Quality', ras.sort_info(ii),savingpath);
% %pngfilename = [num2str(ii,'%02d-'),extractAfter(filenamesavepng,'Cross-Correlation, Auto Correlation and ')];
% filename = ['Spontaneous Activity Analysis for the experiment ',strrep(thisExp.stimulus,'_',' '),', done on ', thisExp.date];
% suptitle_mod(h,filename,2);
% savepngFast(h,savingpath,filename);
% close all;

end