

function Analyze_SpontaneousActivity(datapath, varargin)


totaltime =tic;
if (nargin < 1),    datapath = uigetdir();  end

[thisExp, savingpath] = loadRawData(datapath,'spon','Spontaneous_Activity_Analysis','spondata');

allrasters = CelltoMatUE(thisExp.spiketimes);
binlength = 100/1e3;

clus = thisExp.clusters.goodcells;
endofexp = ceil(nanmax(allrasters,[],'all')/10)*10;

binVec = linspace(0,endofexp,1+endofexp/binlength);
psth = mean(histc(allrasters',binVec),2)*(1/binlength); %#ok

chidx = find (diff(clus(:,1))>0)+1;
%%
h = figure('pos',[200 200 1920 1080],'color','w','vis','on');
subplot_tight(2,1,1,0.05)
rasterPlotter(allrasters',[],'k');
if size(clus,1) > 100
    hold on;
    plot([0 endofexp],[chidx chidx],'r--');
end
ax = gca;
ax.XColor = 'none';
ax.TickDir = 'out';     ax.TickLength = [0.005 0.05];
axis([0 endofexp 0 size(clus,1)+1]);
%yticks(0 : size(clus,1)/6 : size(clus,1))
yticks(chidx(1:floor(length(chidx)/6):end));
yticklabels (clus(chidx(1:floor(length(chidx)/6):end),1));
ylabel('number of electrodes');
title('rasters');

subplot_tight(2,1,2,0.05)
plot(binVec/60,psth,'color',rgb('dodgerblue'),'linewidth',1);
ax = gca;       box off;
ax.TickDir = 'out';     ax.TickLength = [0.005 0.05];
yAx = ceil(max(psth)/2)*2;
axis([0 endofexp/60 0 yAx]);

xticks(0:1:100);
xlabel('time (min)');
yticks(0:yAx/2:yAx);
ylabel('rate (Hz)');
title(['psth (',num2str(binlength*1e3),' ms bins)']);

%[filename, filenamesavepng] = generateRGCname('Cross-Correlation, Auto Correlation and Sorting Quality', ras.sort_info(ii),savingpath);
%pngfilename = [num2str(ii,'%02d-'),extractAfter(filenamesavepng,'Cross-Correlation, Auto Correlation and ')];
filename = ['Spontaneous Activity Analysis for the experiment ',strrep(thisExp.stimulus,'_',' '),', done on ', thisExp.date];
suptitle_mod(h,filename,2);
savepngFast(h,savingpath,filename);
close all;

end