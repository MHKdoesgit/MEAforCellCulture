

function [electstim, mcspulses, h5path , customelecfolder, samplingrate] = electricStimulusTimeStamps(datapath, stimulusname, varargin)

fprintf('reading electrical stimulation events markers....');

mcsh5path = [fileparts(fileparts(which('loadExperimentData.m'))),filesep,'KilosortAddons',filesep,'McsMatlabDataTools'];
if exist(mcsh5path,'dir') 
    if not(exist('McsHDF5.McsData','class'))
        addpath(genpath (mcsh5path));
        import McsHDF5.*
    end
else
    error('There aint no McsMatlabDataTools Toolbox here! how do you want to read the h5 file!, cabr?n!');
end

h5path = [datapath, filesep, 'h5'];
if not(exist(h5path,'dir'))
    h5path = uigetdir(datapath,'Select the Folder with extracted h5 files');
end

h5file = fullfile(h5path,[stimulusname,'.h5']);

cfg = []; cfg.dataType = 'single';
rawrecordeddata = McsHDF5.McsData(h5file,cfg);
samplingrate = round(rawrecordeddata.Recording{1}.AnalogStream{1}.getSamplingRate);

if isempty(rawrecordeddata.Recording{1}.EventStream)
   mcspulses = [];      electstim = [];     customelecfolder = [];
   %electstimh5dat = rawrecordeddata.Recording{1}.EventStream;
   fprintf('no data found\n');
   if ~contains(stimulusname,{'spontaneous','spon'})
   warning('this stimulus aint got any electrical stimulation in it, are you messing around? come back when you know your shit!');
   end
   return;
end

% first load the pre-analyzed electric stim
customelecfolder = [h5path,filesep,'electrical_stimuli'];
customelecfile = [customelecfolder,filesep,stimulusname,'.mat'];
if exist(customelecfile,'file')
    electstim = load(customelecfile);
else
   error('No extracted electrical signal found!, How do you gonna analyze this shit!'); 
end

electstimh5dat = rawrecordeddata.Recording{1}.EventStream{1};

% this is to test and read the raw data
% rawdat = rawrecordeddata.Recording{1}.AnalogStream{1};
% chdat = rawdat.ChannelData;
% chtvec = rawdat.ChannelDataTimeStamps;

electstimevents = electstimh5dat.Events (~cellfun('isempty',electstimh5dat.Events));
electstimlabels = extractAfter (electstimh5dat.Info.Label (~cellfun('isempty',electstimh5dat.Events)),'STG ');
electstimsourcechannel = electstimh5dat.Info.SourceChannelIDs (~cellfun('isempty',electstimh5dat.Events));

for ii = 1: numel(electstimevents)
    
    lb =  lower (strrep (electstimlabels{ii}(3:end),' ','_')); % make the labels for each pulse from the Info structure
    mcspulses.(lb) = double (electstimevents{ii}(1,:)) ./ 1e6; % the time stamps are originally in ?s unit
end
mcspulses.pluseduration = mean(mcspulses.single_pulse_stop - mcspulses.single_pulse_start);
mcspulses.labels = electstimlabels;
mcspulses.sourceChannel = cellfun(@str2double,electstimsourcechannel)';
mcspulses.stimulator = extractBefore(electstimh5dat.Label,';Stimulator');

% toVoltsFac = 10^double(electstimh5dat.Info.Exponent(1));
% mcddat = electstimh5dat.ChannelData(1,:)*toVoltsFac;
fprintf('done\n');
end