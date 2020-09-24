

function Analyze_ElectricalStepStimulus(datapath)

totaltime =tic;
if (nargin < 1),    datapath = uigetdir();  end

thisExp = loadRawData(datapath,'pulse_stimulus','Electrical_Stimulus_Analysis');

if contains(thisExp.stimulus,'sine')
    stimnameending = 'Electrical_stimulus_analysis_for_sine_stimulus';
else
    stimnameending = 'Electrical_stimulus_analysis_for_pulse_stimulus';
end
savingpath = [datapath, filesep,'Data Analysis',filesep,extractBefore(thisExp.stimulus,'_'),'-',stimnameending];

if not(exist(savingpath,'dir'))
    mkdir(savingpath);
end



ft = thisExp.electricalstimulus.single_pulse_start;
%ftoffset = thisExp.electricalstimulus.single_pulse_stop;

out = getRastersandPSTH(ft,thisExp.spiketimes);

for ii = 1: numel(thisExp.spiketimes)
subplot(2,1,1)
cla;
rasterPlotter(out.rastersbefore{ii}',[],'k');
title(ii);
subplot(2,1,2)
cla;
plot(out.psthtime,out.psthbefore(ii,:),'LineWidth',1);

pause;
end




end


%--------------------------------------------------------------------------------------------------%
%----------                               sub-functions                                  ----------%
%--------------------------------------------------------------------------------------------------%

function out = getRastersandPSTH(ft,spk, varargin)

% stim = ft(1:steps:numel(ft));
binLength = 20/1e3;
stimround = round(mean(diff(ft))); % one --|'|--|_| is 5000 ms, 4 frame change
binVec = linspace(0,stimround,1+stimround/binLength);

tbefore = 2;
[rasters, rastersbefore] = deal (cell(numel(spk),1));
[psth, psthbefore] = deal(zeros(numel(spk), length(binVec)));

for ii = 1: numel(spk)
    
    [rasTrials, rasbefore] = deal (cell((length(ft)-1),1));
    eachcellspk = spk{ii};
    
    for jj = 1:numel(ft)-1 %first trial is background intensity
        rasTrials{jj} = eachcellspk(and(eachcellspk > ft(jj),eachcellspk <= ft(jj+1))) - ft(jj);
        rasbefore{jj} = eachcellspk(and(eachcellspk > ft(jj)-tbefore,eachcellspk <= ft(jj+1)-tbefore)) - (ft(jj)-tbefore);
        if isempty(rasTrials{jj}); rasTrials{jj} = NaN; end
    end
    rasters{ii} = CelltoMatUE(rasTrials);
    rastersbefore{ii} = CelltoMatUE(rasbefore);
    
    if size(rasters{ii},2) < 2
        psth(ii,:) = (histc(rasters{ii}',binVec)/size(rasters{ii},1))*(1/binLength);  %#ok to fix situation with 1 spike
    else
        psth(ii,:) = mean(histc(rasters{ii}',binVec),2)*(1/binLength);   %#ok must transpose contspk for histc to work properly.
    end
    
    
    if size(rastersbefore{ii},2) < 2
        psthbefore(ii,:) = (histc(rastersbefore{ii}',binVec)/size(rastersbefore{ii},1))*(1/binLength);  %#ok to fix situation with 1 spike
    else
        psthbefore(ii,:) = mean(histc(rastersbefore{ii}',binVec),2)*(1/binLength);   %#ok must transpose contspk for histc to work properly.
    end
    
end

out.rasters = rasters;
out.rastersbefore = rastersbefore;
out.psth = psth(:,1:end-1);
out.psthbefore = psthbefore(:,1:end-1);
out.psthtime = binVec(1:end-1);
out.binlength = binLength;
out.timebeforespk = tbefore;

end