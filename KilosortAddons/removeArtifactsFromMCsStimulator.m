

function [chdat, pulsestartstop, electpulsedat] = removeArtifactsFromMCsStimulator(eventstream, dat, chanMap, varargin)


fprintf('extracting electrical pulses from the raw data...');

chdat       =    dat.ChannelData(chanMap + 1,:);
chtvec      =    dat.ChannelDataTimeStamps;

activstimulator     = contains (eventstream.Info.Label,'Pulse Start') & ~cellfun('isempty',eventstream.Events)';
activstimulatorstop = contains (eventstream.Info.Label,'Pulse Stop' ) & ~cellfun('isempty',eventstream.Events)';

pulsestart      =   eventstream.Events{activstimulator}(1,:);
pulsestops      =   eventstream.Events{activstimulatorstop}(1,:);
maxdattime      =   max(chtvec);
mindattime      =   min(chtvec);
pulsestart      =   pulsestart(pulsestart > mindattime & pulsestart <= maxdattime);
pulsestops      =   pulsestops(pulsestops > mindattime & pulsestops <= maxdattime);

fonsets         = 1;
madTh           =   40; % start with 40 std and then go to lower one
while numel(fonsets)~= numel(pulsestart)
    madTh           =   madTh-1;
    sumsig          =   abs(sum(chdat,1)); % take absolute sum of all signals
    Thres           =   madTh * mad(diff(sumsig));
    fonsets         =   find( diff(diff(sumsig) >  Thres) <0 ) + 1;
end
stimdur         =   (mean(diff(pulsestart)) / 1e6)- 0.1; % difference between pulses in sec
chgpts          =   diff(fonsets) > (stimdur * dat.getSamplingRate);
stimstrtpts     =   chtvec(fonsets (find(chgpts == 1 ) +1) - 10);% - dat.getSamplingRate/100;
stimendpts      =   chtvec(fonsets (find(chgpts == 1 ) +1) + 5);

if all(stimstrtpts > pulsestart(1:length(stimstrtpts)))
    stimstrtpts     =   chtvec(fonsets ([1, (find(chgpts == 1 ) +1)]) - 10);%  for chunks bigger than 1
    stimendpts      =   chtvec(fonsets ([1, (find(chgpts == 1 ) +1)]) + 5);
end
% choose the best start and end point between the two methods
pulsestart      =   min([pulsestart; stimstrtpts]);
pulsestops      =   min([stimendpts; pulsestops]);

step            =   200;    % check 200 point before the start of the stimulus to get a baseline
stepintime      =   mean(diff(chtvec));     % go forward in time-unit steps
stepsize        =   stepintime;     % amount to go forward in each iteration

% pre-allocation memory for output
pulsestartstop      =   zeros(length(pulsestart),2);
electpulsedat       =   cell(length(pulsestart),1);

tic;
for ii = 1:length(pulsestart)
    
    tbefpulse       =   find (chtvec <= pulsestart(ii),step,'last');
    taftpulse       =   find (chtvec > pulsestops(ii),step,'first');
    
    madbefpulse     =   mad(double(chdat(:,tbefpulse)),0,2); % 0 for mean absolute deviation
    madaftpulse     =   mad(double(chdat(:,taftpulse)),0,2);
    
    chwithpulse     =   (madaftpulse-madbefpulse) > mean(madbefpulse);
    % only include channels that have the pulse
    madbefpulse     =    mad(double(chdat(chwithpulse,tbefpulse)),0,2); % 0 for mean absolute deviation
    madaftpulse     =    mad(double(chdat(chwithpulse,taftpulse)),0,2);
    
    iter            =    0;
    stepintime      =    stepsize;
    while any(abs(madaftpulse-madbefpulse) > abs(mean(madbefpulse))) % cut until the mean of mad
        taftpulse       =   find (chtvec > (pulsestops(ii) + stepintime),stepsize,'first');
        madaftpulse     =   mad(double(chdat(chwithpulse,taftpulse)),0,2);
        stepintime      =   stepintime + stepsize ;
        iter = iter+1;
        if iter > 500, break; end  % break it in case it cannot find the right mad ratio, find better solution to this
    end
    
    tptorem                 =   tbefpulse(end):taftpulse(1);
    pulsestartstop(ii,:)    =   [pulsestart(ii), chtvec(taftpulse(1))];
    electpulsedat{ii}       =   chdat(:,  tptorem);
    
    chdat(chwithpulse, tptorem) = repmat(median(chdat(chwithpulse,tbefpulse),2),1,size(tptorem,2));
    
end
% set ti int16 to match kilosort
chdat = int16(chdat);
fprintf('Took %2.2f seconds...\n',toc);


p12 = (chtvec > pulsestart(ii)-(1.2*tbefpulse(1)) & chtvec < pulsestart(ii+1));

close all
plot(dat.ChannelDataTimeStamps(p12),chdat(12,p12),pulsestart(ii),300,'ro')
hold on
plot(dat.ChannelDataTimeStamps(tbefpulse),chdat(12,tbefpulse),'k','LineWidth',3)
plot(dat.ChannelDataTimeStamps(taftpulse),chdat(12,taftpulse),'m','LineWidth',3)
%n = 7; plot(chtvec,chdat(n,:),'-',stimstrtpts,0,'rx',stimendpts,0,'ko');

end

