

function compareKSsortedtoIgorSorted(ksrawdatapath, igorrawdatapath)


expsks = dir([ksrawdatapath,'/*.mat']);
expsigor = dir([igorrawdatapath,'/*.mat']);

if ~eq(size(expsks,1),size(expsigor,1))
    error('Aint no matching raw data found!');
end
tic;
mintoxc = 20;
xc = cell(size(expsks,1),1);
%prog = Neuropixel.Utils.ProgressBar(15, 'Loading KiloSort dataset: ');
for ii = 1:size(expsks,1)
    
    ksdata = load([expsks(ii).folder,'/',expsks(ii).name]);
    igordata = load([expsigor(ii).folder,'/',expsigor(ii).name]);
    
    xc{ii} = xcorrperexp(ksdata,igordata, mintoxc);
    
end
toc
xcall = cell2mat(xc);

for ii = 1:size(expsks,1)
    subplot_tight(4,4,ii)
    imagesc(xcall(ii).igxcpeak); 
    %title(expsks(ii).name);   
    title(ii)
    set(gca,'ydir','normal');
    axis equal tight off;
end
colormap(gray);

%%
for ii = 1:size(expsks,1)
    subplot(4,4,ii)
    mm = xcall(ii).bestmatch; 
    mm(xcall(ii).peaksim<0.1) = NaN;
    pp = xcall(ii).peaksim; 
    %title(expsks(ii).name);  
    scatter(mm,1:length(mm),pp*40,pp,'filled');
    title(strrep(expsks(ii).name(1:30),'_','-'));
    caxis([0 1]);
    axxy = ceil([size(xcall(ii).ksclus,1) size(xcall(ii).igorclus,1)]/2)*2;
    axis([-1 axxy(1) -1 axxy(2)]);
    xticks(round(0:axxy(1)/4:axxy(1)));
    yticks(round(0:axxy(2)/4:axxy(2)));
end
colormap(flipud(cbrewer('div', 'RdBu',250)));

%%
for ii = 1:size(expsks,1)
    subplot(4,4,ii)
    mm = xcall(ii).bestmatch; 
    mm(xcall(ii).peaksim<0.1) = NaN;
    pp = xcall(ii).igorclus(:,3); 
    %title(expsks(ii).name);  
    scatter(mm,1:length(mm),pp*20,pp,'filled');
    title(strrep(expsks(ii).name(1:30),'_','-'));
    caxis([min(pp) max(pp)]);
    axxy = ceil([size(xcall(ii).ksclus,1) size(xcall(ii).igorclus,1)]/2)*2;
    axis([-1 axxy(1) -1 axxy(2)]);
    xticks(round(0:axxy(1)/4:axxy(1)));
    yticks(round(0:axxy(2)/4:axxy(2)));
end

colormap( [1 0 0; rgb('gold'); rgb('green 3')]);

%%

for ii = 1:size(expsks,1)
    subplot(4,4,ii)
    
    c = [xcall(ii).ksclus(xcall(ii).bestmatch,1:2),nan(size(xcall(ii).igorclus,1),1),xcall(ii).igorclus(:,1:2)];
    c(xcall(ii).peaksim<0.1,:) = NaN;
    
    %mm = xcall(ii).bestmatch; 
    %mm(xcall(ii).peaksim<0.1) = NaN;
    %pp = xcall(ii).igorclus(:,3); 
    pp = xcall(ii).peaksim; 
    %title(expsks(ii).name);  
    scatter(c(:,4),c(:,1),pp*30,pp,'filled');
    title(strrep(expsks(ii).name(1:30),'_','-'));
     caxis([0 1]);
    axxy = [254 254];%ceil([size(xcall(ii).ksclus,1) size(xcall(ii).igorclus,1)]/2)*2;
    axis([-1 axxy(1) -1 axxy(2)]);
    axis square;
    xticks(round(0:axxy(1)/4:axxy(1)));
    yticks(round(0:axxy(2)/4:axxy(2)));
end
colormap(cool);


%%
lagall = [transpose(1:size(igorxc,2))+repmat((lag/max(lag)/3),size(igorxc,2),1),nan(size(igorxc,2),1)]';
for ii =1:size(igorxc,3)

igc = [squeeze(igorxc(:,:,ii));nan(1,size(igorxc,2))];
plot(lagall(:),ii+igc(:))
hold on;
end
axis tight off;
tightfig;



% fpath='E:\Karamanlis_20180712_252MEA10030_mar_sr_le_pc2';
% rstr=load(fullfile(fpath, 'rez.mat'));
% st3=rstr.rez.st3; clear rstr;
% dataIgor=load(['C:\Users\Karamanlis_Dimokrati\Documents\DimosFolder\experiments\'...
%     'Karamanlis_20180712_mar_sr_le_pc2\data_analysis\7_frozencheckerflicker\7_raw_data.mat']);

%%
% maxTime=min([max(igordata.ftimes) 60*10]); %use only 10 minutes of recording
% dt=0.1e-3; %in s
% alltrains = blinkBinner( 0:dt:maxTime,igordata.spiketimes' , 1, 1)'; 
% trainks = blinkBinner( 0:dt:maxTime,ksdata.spiketimes(1)' , 1, 1)'; 
% trainKilosort=gpuArray(trainks);
% 



% 
% %%
% maxLag=6*1e-3/dt;
% idcheck=353;
% indsKilosort=st3(st3(:,2)==idcheck+1 & st3(:,1)<size(alltrains,1),1);
% trainKilosort=zeros(size(alltrains,1),1,'single');
% trainKilosort(indsKilosort)=1;
% trainKilosort=gpuArray(trainKilosort);
% 
% allxcorr=zeros(2*maxLag+1,size(alltrains,2));
% for cellId=1:size(alltrains,2)
%     trainIgor=gpuArray(single(alltrains(:,cellId)));
%     trainxcorr=xcorr(trainKilosort,trainIgor, maxLag,'coeff');
%     allxcorr(:,cellId)=gather(trainxcorr);
% end
% [~,bestMatch]=sort(max(allxcorr),'descend');
% spkIgor=sum(alltrains(:,bestMatch(1)));
% spkKilosort=sum(trainks);%numel(indsKilosort);
% 
% fprintf('Kilosort/Igor spikes : %d/%d (%0.02f), best %d \n',...
%     spkKilosort,sum(spkIgor),spkKilosort/sum(spkIgor),bestMatch(1))
% plot(allxcorr)


end

%--------------------------------------------------------------------------------------------------%
%----------                               sub-functions                                  ----------%
%--------------------------------------------------------------------------------------------------%


function [spkbin, clus] = binexpdata(expData, para)
% loading clusters
clus = expData.clusters;
%bin the spikes for the period of frame timing
%duration_exp=ftimes_ms(end)-ftimes_ms(1);
tvec = 0:para.tbin:para.mstoanalyze; %just take the first 15 minutes, if no it takes too long alternative is; (1)0:tbin:ftimes_ms(end) orm (2)
%the length of the spikes x = 0:tbin:max(sTime_align1);

spkbin = zeros(length(tvec)-1,size(clus,1));
for ii=1:size(clus,1)
    sTime_ms=expData.spiketimes{ii}; %get the spike timing
    %%% align to experiment onset
    sTime_align1=sTime_ms-expData.ftimes(1); %this is done because the zero of x is not the same as for the spikes
    %because the start of the recording doesnt correspond to the start of
    %the frames, so there are spikes before the frames. I want to know how
    %the cross-corr is for when a stimulus is shown not spontanous
    thisbin = histcounts(sTime_align1,tvec);
    spkbin(:,ii)= thisbin(:);
end
end


function xc = xcorrperexp(ksdata,igordata, mintoxc)

mintoanalyze = min([ max(igordata.ftimes)/60, mintoxc]);

%para = expData.stimPara;
para.tbin = 0.5/1e3; %define time bin for histogram
para.mstoanalyze =  mintoanalyze * 60;
para.lagamount = 40/1e3/para.tbin; %we want to go 30 ms
para.lagzoom = 10/1e3/para.tbin;
%para.mea2d = getMEAcoordinatesfromMCD(dp);

[igorbin,xc.igorclus] = binexpdata(igordata, para);
[ksbin,xc.ksclus] = binexpdata(ksdata, para);

[xc.xcmat, xc.lag] = xcorrTwoSpikeTrains(igorbin,ksbin, para.lagamount,1);

xc.igornumspks = sum(igorbin)';


igorxc = zeros(size(xc.xcmat,1),size(xc.xcmat,3),size(xc.xcmat,2));
[peaksim, bestmatch, maxlagloc] = deal(zeros(size(igorbin,2),1));
igxcmat = zeros(size(xc.xcmat,1),size(xc.xcmat,2));
igmaxmat = zeros(size(xc.xcmat,2),size(xc.xcmat,3));

for ii = 1:size(igorbin,2)
   igorxc(:,:,ii) = (squeeze(xc.xcmat(:,ii,:))./xc.igornumspks(ii)); 
   [peaksim(ii) , bestmatch(ii)] = max(max(igorxc(:,:,ii),[],1)); 
   igxcmat(:,ii) =  (igorxc(:,bestmatch(ii),ii));
   
   [~,maxlagloc(ii)] = max(igxcmat(:,ii),[],1);

   igmaxmat(ii,:) = (igorxc(maxlagloc(ii),:,ii));
end

xc.igorxc = igorxc;
xc.igxcbestmatch = igxcmat;
xc.igxcpeak = igmaxmat;
xc.maxlagloc = maxlagloc;
xc.peaksim = peaksim;
xc.bestmatch = bestmatch;

end


