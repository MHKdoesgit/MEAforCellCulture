

function meaFrametimings(varargin)
%
%%% meaFrametimings %%%
%
%
% This function reads the frametimes directly from mcd or msrd file and
% feed them to the new frametimings function. Using this, there is no need
% for data conversion to bin format and the output is tidy.
%
%================================Inputs====================================
%
%   mcdatapath : folder path to the mcd or msrd-h5 file.
%   pulsechannel : channel number that has the pulse data.
%   Show : option to plot figure.
%   Save : option to save the output file as mat.
%   ShapeOffset : pulse shape offset, changes by amplifier.
%   SignalThresh : signal threshold, changes by amplifier.
%   PulseThresh : pulse baseline threshold.
%   MaxHeight : maximum pulse height.
%   Baseline :general baseline value.
%   MSecTick : some bullshit way that mcdata is stored, in 0.1 msec.
%   ZeroADValue2 : baseline of NI AD converter, used in bin conversion only.
%   MicrovoltsPerAD2 : amplified gain to convert to µV.
%   Nblink : value to define onset and off of pulses. Set to 1 if you want
%            to get both onset and offset of each pulse. Otherwise only for stimuli
%            that have the word 1 blink in their name the pulses are extracted
%            whenever the stimuli is changing.
%   SaveFigFormat : option to save the output figure as .fig.
%   OverwriteCriticalParameters : option to overwrite critical values. I
%                                 would not use it if I were you.
%
%================================Output====================================
%
%   frametime : pulses per stimuli change.
%
% written by Mohammad, 07.10.2019, as part of Kilosort functions.

% parse input
p = inputParser();
p.addParameter('mcdatapath', [], @(x) ischar(x));
p.addParameter('pulsechannel', 253, @(x) isnumeric(x));
p.addParameter('Show', true, @(x) islogical(x));
p.addParameter('Save', true, @(x) islogical(x));
p.addParameter('ShapeOffset', 0.1, @(x) isnumeric(x));
p.addParameter('SignalThresh', 40, @(x) isnumeric(x));
p.addParameter('PulseThresh', 0, @(x) isnumeric(x));
p.addParameter('MaxHeight', 500, @(x) isnumeric(x));
p.addParameter('Baseline', 0, @(x) isnumeric(x));
p.addParameter('MSecTick', 0.1, @(x) islogical(x));
p.addParameter('ZeroADValue2', 32768, @(x) isnumeric(x)); % amplifier fix shit, don't change it!
p.addParameter('MicrovoltsPerAD2', 25625/2048, @(x) isnumeric(x)); % amplifier fix shit, don't change it!
p.addParameter('Nblink', [], @(x) isnumeric(x));
p.addParameter('SaveFigFormat', false, @(x) islogical(x));
p.addParameter('OverwriteCriticalParameters', false, @(x) islogical(x));
p.CaseSensitive = false; % make it case-insensitive

p.parse(varargin{:});
params = p.Results;

if isempty(params.mcdatapath) || ~exist(params.mcdatapath,'dir')
    params.mcdatapath = uigetdir([],'Select mc data folder');
end

% get the pulse channel
pulsechannel = params.pulsechannel;

if ~params.OverwriteCriticalParameters
    if pulsechannel> 250 % this starting conditions are very important, otherwise shit will get fucked up!
        % from old rules
        % For 10 KHz
        % metaframes10 ('PulseThresh',10,'SignalThresh',10,'ShapeOffset',0.1,'Invert',true);
        params.PulseThresh = 10;
        params.SignalThresh = 10;
        params.ShapeOffset = 0.1;
        
    elseif pulsechannel < 65  % No invert is need for no bin input!
        % For 25 KHz
        % metaframes25
        % ('PulseThresh',200,'SignalThresh',200,'ShapeOffset',0.04,'Invert',true);
        params.PulseThresh = 200;
        params.SignalThresh = 200;
        params.ShapeOffset = 0.4;
    end
end

% make output folder
ftfolder = [params.mcdatapath,'/frametimes/'];
if ~exist(ftfolder,'dir'), mkdir(ftfolder); end


tt= tic;
% get the data format
[datafilenames, ismcd, ismsrdh5] = mcdormsrdh5(params.mcdatapath);

% set up all the parameters for the frametimings function
ptorm = {'mcdatapath','OverwriteCriticalParameters','pulsechannel'};
frops = cell(1,(numel(p.Parameters)-numel(ptorm))*2);
frops(1:2:end) = p.Parameters(~contains(p.Parameters, ptorm));
frops(2:2:end) = struct2cell(rmfield(params,ptorm));

for ii = 1: numel(datafilenames)
    
    mcname = datafilenames{ii};
    
    if ismcd
        mcdat = readMCDanalongData(params.mcdatapath, mcname,pulsechannel);
    end
    
    % get the blinks right
    frops = setnblinks(mcname, frops, params);
    
    % here is essentially is the old metaframe function
    fprintf('Extracting frame timings from recording ''%s''...\n',mcname);
    % call frametimings
    frametimings('DirectInput', mcdat,'FileName',mcname,frops{:});
    % move and rename files because of bullshit dairy function inside frametimings
    fprintf('Moving files to ''%s''...\n',ftfolder);
    
    fileext = {'log','mat','fig','png'};
    for ex=1:length(fileext)
        srcfile =[mcname(1:end-4),'_frametimings.',fileext{ex}];
        if exist(srcfile,'file')
            movefile(srcfile,[ftfolder,mcname(1:end-4),'_frametimings.',fileext{ex}]);
        end
    end
    disp('Done.');
end

if ismcd
    clear mexprog; %unload DLL
end
toc(tt);

end

%--------------------------------------------------------------------------------------------------%
%----------                               sub-functions                                  ----------%
%--------------------------------------------------------------------------------------------------%

function [outputfilename, ismcd, ismsrdh5, mcddll] = mcdormsrdh5(dp)

mcdfilenames = dir([dp,filesep,'*.mcd']);
h5filenames = dir([dp,filesep,'*.h5']);

ismcd = ~isempty(mcdfilenames);
ismsrdh5 = ~isempty(h5filenames);

if ~ismcd && ~ismsrdh5
    error('There aint no valid data shit in this folder, come at me with mcd or msrd-h5 format!');
end

if ismcd && ismsrdh5
    error('DAFAQQ, you have recording from both amplifiers in this folder, Shit is confusing');
end

if ismsrdh5
    [~, reindex]=sort(str2double(regexp(({h5filenames(:).name}),'\d+','match','once')));
    h5filenames={h5filenames(reindex).name}';
    outputfilename = h5filenames;
    mcddll = [];
end

if ismcd
    [~, reindex]=sort(str2double(regexp(({mcdfilenames(:).name}),'\d+','match','once')));
    mcdfilenames={mcdfilenames(reindex).name}';
    outputfilename = mcdfilenames;
    %load the dll file
    [dllpath,libtoload] = getMCSdllPath();
    mcddll=mexprog(18, [dllpath, filesep, libtoload]);  %set dll library
end

end

%--------------------------------------------------------------------------------------------------%

function [mcddat, stimsamples] = readMCDanalongData(dp, mcname, pulsechannel)

mcdfile = [dp,filesep,mcname]; %get mcd path
[~, hfile] = mexprog(1, mcdfile); %open file
[~, mcdfileInfo] = mexprog(3, hfile); %get file info
% in case of loading all analog channels comment out these line

% [~, chinfos] = mexprog(4, hfile,0:(mcdfileInfo.EntityCount-1)); %get channel info
% labellist = {chinfos.EntityLabel};  %extract labels of the entities
% [~, analogchans] = getChannelMapForRawBinary(labellist);
% [~, volinfos] = mexprog(7, hfile,analogchans); % voltage range for the pulse channels are different

[~, volinfos] = mexprog(7, hfile, pulsechannel); % voltage range for the pulse channels are different
maxVoltage=volinfos(1).MaxVal; minVoltage=volinfos(1).MinVal;
%resVoltage=volinfos(1).Resolution;
newRange = 2^15*[-1 1];
multFact=range(newRange)/(maxVoltage-minVoltage);

stimsamples = mcdfileInfo.TimeSpan/mcdfileInfo.TimeStampResolution;
[~,~,mcddat]=mexprog(8,hfile, pulsechannel-1,0,  stimsamples+1);%read data
if sum(mcddat)==0 % dumm trick to guess the correct length of the data, in case last index is missing
    [~,~,mcddat]=mexprog(8,hfile, pulsechannel-1,0,  stimsamples);%read data
end
% re-scale the output
mcddat = (mcddat.*multFact)';
%close mcd file
mexprog(14, hfile);

end

%--------------------------------------------------------------------------------------------------%

function [dllpath,libtoload] = getMCSdllPath()
%GETMCSDLLPATH Summary of this function goes here
dlllocation = which('load_multichannel_systems_mcd');
dllpath = fileparts(dlllocation);

switch computer()
    case 'PCWIN'; libtoload = 'nsMCDLibraryWin32.dll';
    case 'GLNX86'; libtoload = 'nsMCDLibraryLinux32.so';
    case 'PCWIN64'; libtoload = 'nsMCDLibraryWin64.dll';
    case 'GLNXA64'; libtoload = 'nsMCDLibraryLinux64.so';
    case 'MACI64'; libtoload = 'nsMCDLibraryMacIntel.dylib';
    otherwise
        disp('Your architecture is not supported'); return;
end
end

%--------------------------------------------------------------------------------------------------%

function ops = setnblinks(mcname, ops, params)
% here is essentially is the old metaframe function
Nblink = str2double(mcname(1,regexp(lower(mcname),'blink')-1));
if isempty (Nblink)
    Nblink =str2double(mcname(1,cell2mat(regexp(lower(mcname),'blinck'))-1));
end
% the index of nblink cell from the ops cell
inputNblinkidx = find(contains(lower(ops(1:2:end)),'nblink'))*2;

if isempty (Nblink) || isnan (Nblink)
    if isempty(ops{inputNblinkidx})
        %disp('No Nblinks, You Shall have Nblinks of 2');
        Nblink = 2;  % Default value
    else
        Nblink = params.Nblink;
    end
end
ops{inputNblinkidx}= Nblink;

end
