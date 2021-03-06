

function spikesfromElectricStimulation(datapath)

fprintf('reading electrical stimulation events markers....');
%%
mcsh5path = [fileparts(which('loadExperimentData.m')),filesep,'McsMatlabDataTools'];
if exist(mcsh5path,'dir')
    addpath(genpath (mcsh5path));
    import McsHDF5.*
else
    error('There aint no McsMatlabDataTools Toolbox here! how do you want to read the h5 file!, cabr?n!');
end

h5path = [datapath, filesep, 'h5'];
if not(exist(h5path,'dir'))
    h5path = uigetdir(datapath,'Select the Folder with extracted h5 files');
end

h5file = fullfile(h5path,[stimulusname,'.h5']);

cfg = []; cfg.dataType = 'raw';
rawrecordeddata = McsHDF5.McsData(h5file,cfg);
samplingrate = round(rawrecordeddata.Recording{1}.AnalogStream{1}.getSamplingRate);

if isempty(rawrecordeddata.Recording{1}.EventStream)
   electstim = []; 
   electstimh5dat = rawrecordeddata.Recording{1}.EventStream;
   fprintf('no data found\n');
   if ~contains(stimulusname,{'spontaneous','spon'})
   warning('this stimulus aint got any electrical stimulation in it, are you messing around? come back when you know your shit!');
   end
   return;
end
electstimh5dat = rawrecordeddata.Recording{1}.EventStream{1};

% this is to test and read the raw data
rawdat = rawrecordeddata.Recording{1}.AnalogStream{1};
chdat = rawdat.ChannelData;
chtvec = rawdat.ChannelDataTimeStamps;


electstimevents = electstimh5dat.Events (~cellfun('isempty',electstimh5dat.Events));
electstimlabels = extractAfter (electstimh5dat.Info.Label (~cellfun('isempty',electstimh5dat.Events)),'STG ');
electstimsourcechannel = electstimh5dat.Info.SourceChannelIDs (~cellfun('isempty',electstimh5dat.Events));

for ii = 1: numel(electstimevents)
    
    lb =  lower (strrep (electstimlabels{ii}(3:end),' ','_')); % make the labels for each pulse from the Info structure
    electstim.(lb) = double (electstimevents{ii}(1,:)) ./ 1e6; % the time stamps are originally in ?s unit
end
electstim.pluseduration = mean(electstim.single_pulse_stop - electstim.single_pulse_start);
electstim.labels = electstimlabels;
electstim.sourceChannel = cellfun(@str2double,electstimsourcechannel)';
electstim.stimulator = extractBefore(electstimh5dat.Label,';Stimulator');

elstimtvec = electstim.single_pulse_start .* 1e6;
pulsespkdur = 3;
pidxlen = ((pulsespkdur*samplingrate)/1e2)+1;
pulses = false(length(elstimtvec), size(chtvec,2));

for jj = 1: numel(electstim.single_pulse_start)
    pulsereg = find (double(chtvec) >= elstimtvec(jj) & double(chtvec) <= (elstimtvec(jj) + pulsespkdur*1e4));
    if length(pulsereg) > pidxlen
        pulsereg = pulsereg(1:pidxlen);
    elseif length(pulsereg) < pidxlen
        pulsereg = pulsereg(1):pulsereg(1)+(pidxlen-1);
    end
    pulses(jj, pulsereg) = true;   
end
chdatpulses = chdat(:,any(pulses,1));
chtvecpulses = double(chtvec(any(pulses,1))) ./ 1e6;
%%
figure
for ii = 1:60
pulsespks = reshape(chdatpulses(ii,:),[],numel(electstim.single_pulse_start))';
cla;
subplot(1,2,1)
plot(pulsespks');
xlim([0 100]);
ax = gca;
title(ii);
subplot(1,2,2)
plot(diff(pulsespks)')
axis([ax.XLim, ax.YLim]);
pause;
end
%%
plot(double(chtvec)/1e6,chdat(14,:));
hold on;
plot(electstim.single_pulse_start,max(chdat(12,:))-range(chdat(12,:))/8,'ro')
plot(electstim.single_pulse_stop,max(chdat(12,:))-range(chdat(12,:))/8,'kx')
plot(chtvecpulses, chdatpulses(12,:),'k.')

end