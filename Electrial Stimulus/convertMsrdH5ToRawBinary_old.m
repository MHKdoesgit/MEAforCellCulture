

function [bininfo] = convertMsrdH5ToRawBinary(ops, targetpath)

% status = system('"Multi Channel DataManager.exe" &')
% status = system('set PATH=' 'C:\Program Files\ Multi Channel DataManager' ' && ''Multi Channel DataManager.exe &');
% status = system('set path=%path:"C:\Program Files\Multi Channel DataManager\";=% & "Multi Channel DataManager.exe" &');
% dataHFfilename = '2019-05-21T18-13-48HippocampalCultures_ControlReal_spontaneous.h5';
% cfg = [];
% % cfg.dataType = 'single';
% cfg.dataType = 'raw';
% data = McsHDF5.McsData([datapath,'/',dataHFfilename],cfg);
% %%
% % One can convert data loaded with the 'raw' option to meaningful units
% % either manually (in this example for the first channel of an analog stream):
% %
% % converted_data = (data.Recording{1}.AnalogStream{1}.ChannelData(1,:) - ...
% %     double(data.Recording{1}.AnalogStream{1}.Info.ADZero(1))) * ...
% %     double(data.Recording{1}.AnalogStream{1}.Info.ConversionFactor(1));
%
%
% raw_data = data.Recording{1}.AnalogStream{1}.getConvertedData(cfg);
% filtered_data = data.Recording{1}.AnalogStream{2}.getConvertedData(cfg);
% samplerate = data.Recording{1}.AnalogStream{2}.getSamplingRate;

%--------------------------------------------------------------------------
tic;
% get msrd filenames
h5filenames         =   dir([ops.root,filesep,'*.h5']);
[~, reindex]        =   sort(str2double(regexp(({h5filenames(:).name}),'\d+','match','once')));
h5filenames         =   {h5filenames(reindex).name}';
Nfiles              =   numel(h5filenames);
%--------------------------------------------------------------------------
% get information about the recording time
% config files for dataloading
cfg                 =   [];
cfg.dataType        =   'raw';

stimdata            =   cell(numel(h5filenames),1);
stimsamples         =   zeros(numel(h5filenames),1);
elecstimflag        =   zeros(numel(h5filenames),1);

for imcd    =   1 : numel(h5filenames)
    
    h5pathname      =   [ops.root,filesep,h5filenames{imcd}]; %get mcd path
    stimdata{imcd}  =   McsHDF5.McsData(h5pathname,cfg);
    stimsamples(imcd)   =   size(stimdata{imcd}.Recording{1}.AnalogStream{1}.ChannelDataTimeStamps,2);
    elecstimflag(imcd)  =   ~isempty(stimdata{imcd}.Recording{1}.EventStream);
    
end
bininfo.stimsamples     =   stimsamples;

streamtype          =    cell(size(stimdata{imcd}.Recording{1}.AnalogStream));
for ii  =    1 : size(stimdata{imcd}.Recording{1}.AnalogStream,2)
    streamtype{ii}  =    stimdata{imcd}.Recording{1}.AnalogStream{ii}.Label;
end
filteredstream      = (contains(streamtype,'Filter')); % get only the filtered stream

H5fileInfo          =   stimdata{imcd}.Recording{1}.AnalogStream{filteredstream}.Info;
NchanTOT            =   size(H5fileInfo.ChannelID,1);
bininfo.NchanTOT    =   NchanTOT;
fs          =    stimdata{imcd}.Recording{1}.AnalogStream{filteredstream}.getSamplingRate;
bininfo.fs  =    fs; % sampling frequency
fprintf('Total length of recording is %2.2f min...\n',sum(stimsamples)/fs/60);
%--------------------------------------------------------------------------
% get the channel names based on the map of the array
labellist           =   {H5fileInfo.Label};
chanMap             =   getChannelMapForRawBinary(labellist,'dataformat','msrd','channelnumber',NchanTOT);
%--------------------------------------------------------------------------
fprintf('Saving .mcd data as .dat...\n');

maxSamples          =    120e5;% chunk size

fidOut      =    fopen(targetpath, 'W'); %using W (capital), makes writing ~4x faster
if  any(elecstimflag)
    savingpath      =    [ops.root,filesep,'electrical_stimuli'];
    if not(exist(savingpath, 'dir')); mkdir(savingpath); end
end

for iFile   =   1 : Nfiles
    
    h5dat       =    stimdata{iFile}.Recording{1}.AnalogStream{filteredstream};
    nsamples    =    stimsamples(iFile);
    Nchunk      =    ceil(nsamples/maxSamples);
    %a = cell(1,Nchunk);
    
    if elecstimflag(iFile)
        elecstimdat     =   stimdata{iFile}.Recording{1}.EventStream{1};
        %[~,nelectevents]=   (cellfun(@size,elecstimdat.Events));
        %nelectevents    =   max(nelectevents);
        electpulsetimes =   cell(NchanTOT,Nchunk);
        electpulsedat   =   cell(NchanTOT,Nchunk);
        %electpulsetimes =   zeros(nelectevents,2);
        %electpulsedat   =   cell(nelectevents,1);
        %elecidx         =   0;
        
        % if we have only one pulse and the delay between the start pulse and stop pulse is more than
        % one minute, then something is suspicious. This is all because the stupid amplifier doesn't
        % make pulseses with custom made stimuli.
        
        %if size(elecstimdat.Events{1},2)==1 && (elecstimdat.Events{2}(1,1) - elecstimdat.Events{1}(1,1) > 60e6 )
        %   [elecstimdat, stimproperties ] = pulsesforCustomElectricStimuli(ops.root, h5filenames{iFile}, elecstimdat);
        %   electrical_stimuli.stimproperties = stimproperties;
        %end
    end
    
    for iChunk  =   1 : Nchunk
        offset          =   max(0, (maxSamples * (iChunk-1)));
        sampstoload     =   min(nsamples-offset,maxSamples);
        lastidxtoload   =   offset+sampstoload; % added by MHK to avoid end index crashes
        if lastidxtoload > size(h5dat.ChannelDataTimeStamps,2)
            lastidxtoload   =    size(h5dat.ChannelDataTimeStamps,2);
        end
        cfg.window  =   double([h5dat.ChannelDataTimeStamps(offset+1) h5dat.ChannelDataTimeStamps(lastidxtoload)]) / 1e6;
        % to convert from microseconds to sec for more info check McsHDF5.TickToSec
        
        dat         =   h5dat.readPartialChannelData(cfg);
        if elecstimflag(iFile)
            [dat, pulsestartstop, electpulse, electpulsepara]   =    removeElectricalStimArtifacts(elecstimdat, dat, chanMap);
            %electpulsetimes(elecidx+1: elecidx+size(pulsestartstop,1),:)    =    pulsestartstop;
            %electpulsedat(elecidx+1: elecidx+size(pulsestartstop,1))        =    electpulse;
            electpulsetimes(:,iChunk)                           =    pulsestartstop;
            electpulsedat (:, iChunk)                           =    electpulse;
            %elecidx     =   elecidx + size(pulsestartstop,1);
        else
            dat         =   int16(dat.ChannelData(chanMap + 1,:));
        end
        %a{iChunk} = dat;
        fwrite(fidOut, dat, 'int16');
    end
    
    if elecstimflag(iFile)
        % to rearrage all the pulses and thier time together
        [r,~]      =       cellfun(@size,electpulsedat);
        r          =       unique(r,'rows');
        numelectpulses     =   sum(r(1,:));
        [elpldat,elpltime] =   deal(nan(NchanTOT, numelectpulses, electpulsepara.nt0));
        iter       =       0;
        for ii     =   1:Nchunk
            td     =   reshape(cell2mat(electpulsedat(:,ii)),r(1,ii),NchanTOT,[]);
            td     =   permute(td,[2 1 3]);
            tt     =   reshape(cell2mat(electpulsetimes(:,ii)),r(1,ii),NchanTOT,[]);
            tt     =   permute(tt,[2 1 3]);
            elpldat(:,iter+1:iter+r(1,ii),:)       =   td;
            elpltime(:,iter+1:iter+r(1,ii),:)      =   tt;
            iter   =   iter + r(1,ii);
        end
        electrical_stimuli.pulsetimes   =    elpltime./1e6;
        electrical_stimuli.pulsedata    =    elpldat;
        electrical_stimuli.para         =    electpulsepara;
        
        save([savingpath,filesep,extractBefore( h5filenames{iFile},'.h5'),'_electrical_stimuli.mat'],...
            '-v7.3','-struct','electrical_stimuli');
        
        clearvars pulsestartstop electpulse electpulsepara electpulsetimes electpulsedat r numelectpulses...
            elpldat elpltime iter ii td tt;
    end
    
    %report status
    fprintf('Time %3.0f min. Mcd files processed %d/%d \n', toc/60, iFile,Nfiles);
end
% fclose(fidOut);
%--------------------------------------------------------------------------
end

