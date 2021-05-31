

function [ experiment ] = loadExperimentData( datapath , varargin)
%
%%% loadExperimentData %%%
%
%
% This function loads everthing about on expriment folder, this include all
% parameters of stimulus, goodcells and clusters, frametimings, spiketimings
% and binned spikes.
%
% ===============================Inputs====================================
%
%    dataPath : path to experiment folder
%
%================================Output====================================
%
%   exp : multi-dimentional structure containing parameters of stimulus,
%         goodcells and clusters, frametimings, spiketimings
%         and binned spikes.
%
% based on loadExperiment function from Fernando,
% modified with new functions by Mohammad, 30.07.2014
% update to make indiviual exp files per experiment on 14.06.2015
% update the saving function for structures on 09.03.2016
% added new options for loading clusters with lower scores on 02.01.2017.
% added new options to load KiloSort data on 18.06.2019
% added new options for KS2 and phy2 on 19.02.2020.

if nargin < 1
    datapath = uigetdir();
end

tic;
if (any(datapath == 0) || isempty(datapath) || ~isfolder(datapath))
    error('You must supply a dataPath path');
end

bestScore = 1;  % by default only the cells with score 1 are considerd as good.
worstScore = 4;

if nargin == 2  % if only one input is given that is considered as worst score and 1 as best score.
    bestScore = 1;
    worstScore = varargin{1};
end

if nargin == 3    % both best and worst scores are defined by user.
    bestScore = varargin{1};
    worstScore = varargin{2};
end

msrdpath = [datapath,filesep,'msrd'];
h5path = [datapath,filesep,'h5'];
kspath = [h5path,filesep,'ks_sorted'];

if not(exist(msrdpath,'dir')) 
    error('there aint no data here, give the correct path with msrd folder in it'); 
end

if not(exist(h5path,'dir')) 
    error('there aint no data here, give the correct path with h5 folder in it'); 
end

if not(exist(h5path,'dir')) 
    error('there aint no spike sorted data here, you shold have a ks_sorted folder here!'); 
end

% this part is to load correct rasters files
if not(exist([kspath,filesep,'ksrasters',filesep,'ksrasters.mat'],'file')) 
    ksdata = readKilosortData(kspath);
else
    ksdata = load([kspath,filesep,'ksrasters',filesep,'ksrasters.mat']);   
end


expSharedInfo.originalFolder = datapath;
[~, expName] = fileparts(expSharedInfo.originalFolder);
expSharedInfo.expName = expName;
expSharedInfo.date = datemaker(datapath);

% making folder for saving data
if not(exist(fullfile(datapath,'Data Analysis','Raw Data'),'dir'))
    mkdir(fullfile(datapath,'Data Analysis','Raw Data'));
end

fprintf('loading good cells...');

clusidx = ksdata.clusters(:,3) <= worstScore & ksdata.clusters(:,3) >= bestScore;
clustersGoodCells = ksdata.clusters(clusidx,:);
fprintf('done!,took... %.1f sec\n',toc);

[stimulipara, stimulinames] = loadStimulusParameters(datapath);

for ii  = 1: numel(stimulinames)
    
    thisExp.originalfolder          =   expSharedInfo.originalFolder;
    thisExp.experimentname          =   expSharedInfo.expName;
    thisExp.stimulus                =   stimulinames{ii};
    thisExp.clusters.goodcells      =   clustersGoodCells;
    thisExp.clusters.allclusters    =   ksdata.clusters;
    thisExp.clusters.goodIndex      =   clusidx;
    
    thisExp.amplitudes              =   ksdata.amplitudes(clusidx,ii);
    thisExp.spiketimes              =   ksdata.spike_times(clusidx,ii);
    temptemp                        =   ksdata.template_info(clusidx);
    thisExp.templates               =   {temptemp.templates}';
    [thisExp.electricalstimulus, thisExp.pulseinfo, h5path, electinfopath] = electricStimulusTimeStamps(datapath, stimulinames{ii});
    
    if isempty(thisExp.electricalstimulus)
        thisExp = rmfield( thisExp, 'electricalstimulus');
    end
    
    thisExp.mea.numchannels         =   ksdata.sort_params.params_py.n_channels_dat;
    thisExp.mea.map                 =   ksdata.sort_params.channel_positions;
    thisExp.mea.labels              =   ksdata.sort_params.channel_map;
    thisExp.stimPara                =   stimulipara{ii};
    thisExp.sortinginfo             =   ksdata.sort_info(clusidx);
    thisExp.stim_start_end          =   ksdata.sort_params.stim_start_end;
    thisExp.samplingrate            =   ksdata.sort_params.sampling_rate;
    thisExp.h5path                  =   h5path;
    thisExp.electstimpath           =   electinfopath; 
    thisExp.date                    =   expSharedInfo.date;
    
    experiment = thisExp;
    save([datapath,'/Data Analysis/Raw Data/',stimulinames{ii},' for experiment on ',thisExp.date,'.mat'],'-v7.3','-struct','thisExp');
    
    fprintf(['pre-analyzing stimulus ',stimulinames{ii}, '...is finito, took... %.1f sec\n'],toc);
    clearvars thisExp ;
end

end


