

function wnelectstim = whiteNoiseElectricalStimulus(mins, stimfrq, pulsedur, seed, savingpath)

if nargin < 5, savingpath = []; end

msstimdur = mins * 60* 1e3;

numstimchange = floor(msstimdur ./ (1e3/stimfrq));
meanStim = 0;
stdStim = 500;
%pfr = 0;

if seed >= 0, error('What is this shit, only negative seed!'); end
stimraw = meanStim + stdStim * (gasdev(seed, numstimchange));

pulsefreq = (1/stimfrq)*1e6; % amplifier units are in ?S
iter = pulsefreq;       % start one stimfrq late

stimpulse = cell(1,numstimchange);
for ii = 1:numstimchange
    %origgsdvvalue = ((stimraw(ii)-meanStim) / stdStim);
    if stimraw(ii) < 0 %origgsdvvalue < 0
        pulsetype = 'anodic';
    else
        pulsetype = 'cathodic-anodic';
    end
    
    stimpoint = generateElectricPulses(pulsedur,stimraw(ii),0,pulsetype);
    %tdelays = iter - max(stimpoint(1,:));
    tdelays = iter;
    
    stimpulse{ii} = [[tdelays;0],[tdelays+1;0]+stimpoint];
    
    iter = iter+pulsefreq;
    
end

stimpulseall = cell2mat(stimpulse);
stimpulseall(2,:) = stimpulseall(2,:) * 1e3; % convert to milivolt
stimpulseall = int64(stimpulseall);

recdur = double(max(stimpulseall(1,:))/60e6);
fprintf('total required recording time:\t\t%0.1f minutes (%.2f hours).\n',recdur,recdur/60);


wnelectstim.stimulus = double(stimpulseall);
wnelectstim.para.meanstim = meanStim;
wnelectstim.para.stdstim = stdStim;
wnelectstim.para.stimfreq = stimfrq;
wnelectstim.para.numstimchange = numstimchange;
wnelectstim.para.seed = seed;
wnelectstim.para.stimdurationms = msstimdur;


if ~isempty(savingpath)
    wrfilename = sprintf('electricWhiteNoiseStimulus_mean%.1f_std%.f_at%.2fHz',meanStim,stdStim,stimfrq);
    writeMEAstimulusTxt(stimpulseall', savingpath, wrfilename);
    % save everthing to mat file to be sure!
    save([savingpath,filesep,wrfilename,'.mat'],'-v7.3','-struct','wnelectstim');
end


%
% stimrawInt = int64(stimraw);
% stimtimes = 0: (1e3/stimfrq): msstimdur;
% stimtimes = stimtimes(1:numstimchange)';
%
% stimtimesmicrosec = int64(stimtimes * 1e3) + pfr;
% if stimtimesmicrosec(1) == 0, stimtimesmicrosec = stimtimesmicrosec + 1; end % to avoid negtive time
% stimtimebsf = transpose([stimtimesmicrosec-1, stimtimesmicrosec, stimtimesmicrosec+pulsedur+1]);
%
% stimVoltCurr = transpose([zeros(size(stimrawInt)), stimrawInt, zeros(size(stimrawInt))]);
%
% stimtowrite = [stimtimebsf(:),stimVoltCurr(:)];
%


end
