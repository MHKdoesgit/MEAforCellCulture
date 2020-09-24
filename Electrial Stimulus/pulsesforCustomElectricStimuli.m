

function [esdatout, stimout] = pulsesforCustomElectricStimuli(h5path, h5name, esdat)

stimname = [extractBefore(h5name, '.h5'),'.txt'];
if exist([h5path, filesep, stimname],'file')
    stimpath = [h5path, filesep, stimname];
else
    stimpath = [fileparts(h5path),filesep,'stimuli',filesep,stimname];
end

if not(exist(stimpath,'file'))
    error('the fuck! there is no text file for the stimulus here, what do expect me to do!!');
end

ogstimpath = [extractBefore(h5path,[filesep,'MEA',filesep]),[filesep,'MEA',filesep,'Electrical Stimulus']];
if not(exist(ogstimpath,'dir'))
    ogstimpath = stimpath;
end

ogstimfilename = [cell2mat(extractBetween(stimname,'_', '.txt')),'.mat'];
ogstimpath = [ogstimpath,filesep,ogstimfilename];

if not(exist(ogstimpath,'file'))
    error('There aint no stimulus file here, how do expect me to figure out stimulus properties?! cabron!');
end

ogstim = load(ogstimpath);

fid = fopen(stimpath,'r');
es = textscan(fid,'%d64%d64','Delimiter','\t');
fclose(fid);

estim = double(cell2mat(es));

activstimulator     = contains (esdat.Info.Label,'Pulse Start') & ~cellfun('isempty',esdat.Events)';
activstimulatorstop = contains (esdat.Info.Label,'Pulse Stop' ) & ~cellfun('isempty',esdat.Events)';

pulsestart      =   esdat.Events{activstimulator}(1,:);
pulsestops      =   esdat.Events{activstimulatorstop}(1,:);

lastestimepoint = estim(end,1);

if lastestimepoint > pulsestops
    error('somehow, your last pulse comes after the end of electrical stimulation, good new! your recording is good for shit!');
end

stimnamefrq = str2double(cell2mat(extractBetween(stimname,'_at','Hz')));
stimfrq = 1/ ((estim(2,1)-estim(1,1))/1e6);

if ~isequal(stimnamefrq, stimfrq)
    warning('hmmm, something looks fishy, no update frequency is mentioned in the file name');
end

stiminsec = (1/stimfrq) * 1e6; % from micro-seconds to sec


newpulsetime = (rem(estim(:,1),stiminsec)==0); % this is when each new pulse is happening

newpulsestart = estim(newpulsetime,1);
newpulsestart = newpulsestart(2:end); % the fist data point is zero

%ogstim
[ca, cb, cc] = ndgrid(ogstim.pulsedurations, ogstim.pulseamplitudes, ogstim.pulsetypes);
possiblecombs = [ca(:), cb(:), cc(:)];
numcombs = size(possiblecombs,1);

stimrepeat = cell(numcombs,1);
for ii = 1:ogstim.numrepeats
    
    randstimorder =  psudoRandPerm(ogstim.seeds(ii), 1:numcombs);
    stimrepeat{ii} = possiblecombs(randstimorder,:);
    
end
stimall = cell2mat(stimrepeat);
stimdur = stimall(:,1);
ancathstims = (stimall(:,3) < 3 | stimall(:,3) > 4); % 3 and 4 are one-sided pulses
stimdur(ancathstims) = (stimdur(ancathstims) * 2) + 3; % add ones for each step (required by stupid amplifier)
stimdur(~ancathstims) = stimdur(~ancathstims) + 2;

newpulsestart = newpulsestart + 1;
newpulsestop = newpulsestart + stimdur;

pulsesstartdurrecording = int64(newpulsestart) + pulsestart;
pulsesstopdurrecording  = int64(newpulsestop) + pulsestart;


eventsstart = [pulsesstartdurrecording';repmat(esdat.Events{1}(2:end,1),1,size(stimall,1))];
eventsstop = [pulsesstopdurrecording';repmat(esdat.Events{2}(2:end,1),1,size(stimall,1))];

%eventsstart = [esdat.Events{1}, eventsstart];
%eventsstop = [eventsstop, esdat.Events{2}];

esdatout.Events = esdat.Events;
esdatout.TimeStampDataType = esdat.TimeStampDataType;
esdatout.Label = esdat.Label;
esdatout.Info = esdat.Info;


esdatout.Events{1} = eventsstart;
esdatout.Events{2} = eventsstop;

stimout = ogstim;
stimout.stimcombinations = stimall;
stimout.stimdurations = stimdur;
stimout.possiblecombs = possiblecombs;

%
% plot(estim(:,1),estim(:,2))
% hold on
% plot(newpulsestop,eschpt,'kx')
% plot(newpulsestart,eschpt,'r.')


end