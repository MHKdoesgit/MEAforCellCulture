

function [chdat, pulsestartstop, electpulsedat, para] = removeElectricalStimArtifacts(eventstream, dat, chanMap, varargin)

if nargin > 3 
    if strcmpi(varargin{1},'old')
       [chdat, pulsestartstop, electpulsedat] = remelecpulses_old(eventstream, dat, chanMap);
       return;
    end
end

fprintf('extracting electrical pulses from the raw data...');

tic;
chdat       =    dat.ChannelData(chanMap + 1,:);
chtvec      =    dat.ChannelDataTimeStamps;
chrawdat = single(chdat)';

% some values to get correct electrical signals
para.local_range       =   [5  0];
para.nt0               =      300;
para.nt0min            =       75;
para.threshold         =    -1500;
[nT, nChan] = size(chrawdat);

% mean with blanking, from kilosort my_min
smin     =      single(pulseminpoint(chrawdat, para.local_range, [1 2]));
ppeaks   =      single(chrawdat < smin+1e-2 & chrawdat < para.threshold);

% this is to correct error in detect too many close by pulses
for ii = 1: nChan
    p = ppeaks(:,ii);
    ppts = find(p);
    wrngpls = ppts (find(diff(ppts) < 1e3) +1);
    ppeaks(wrngpls,ii) = 0;
end
% this is to unify the number of pulses among all channels
npulseperch = sum(ppeaks,1);
[npl,nplbins] = histcounts(npulseperch,[unique(npulseperch), Inf]);
nplbins = nplbins( nplbins < Inf);
numchlim = 10;
numpulses = max (nplbins(npl > numchlim));  % in case  it missed some pulses in many channels, take the most probable one
if isempty(numpulses)
    numchlim = 5;
    numpulses = max (nplbins(npl > numchlim));    
end
%numpulses = mode(npulseperch);
chwithcorrpulse = find (npulseperch == numpulses);
for ii = 1: nChan
    if npulseperch(ii) < numpulses || npulseperch(ii) > numpulses
        chtocmp = ppeaks ( : , chwithcorrpulse (randperm(numel(chwithcorrpulse),numchlim)));
        ppts = find ( mean(chtocmp,2) > 0.1);
        % replace the location of missing pulses in noisy channels with the good pulses
        corrpulses = [ppts(not(diff(ppts) < 1e3) ); ppts(end)];
        if length(corrpulses) < numpulses
            ppts = find ( mean(chtocmp,2) > 0.01);
            corrpulses = [ppts(not(diff(ppts) < 1e3) ); ppts(end)];
        end  
        p =  zeros(size(ppeaks(:,ii)));
        p(corrpulses) = 1;
        ppeaks(:,ii) = p;
    end
end

[row, col, ~] = find(ppeaks);
%mu = - mu;

% n = 1;plot(1:size(S1,1),S1(:,n),row(col==n),mu(col==n),'ro')
% hold on
% plot(row(1)-25:row(1)+25,S1(row(1)-25:row(1)+25,1),'m','LineWidth',3)
%for n=1:60; plot(1:size(chrawdat,1),chrawdat(:,n),row(col==n),mu(col==n),'ro'); title(n); pause; end


% times around the peak to consider
dt = 1:para.nt0;
dc = 0;

% temporal indices
indsT = repmat((row-para.nt0min)', numel(dt), 1) + repmat(dt', 1, numel(row));
% spatial indices (here is zero)
indsC = repmat(col', numel(dc), 1) + repmat(dc', 1, numel(col));

indsC(indsC<1)     = 1;
indsC(indsC>nChan) = nChan;

indsT = permute(indsT, [1 3 2]);
indsC = permute(indsC, [3 1 2]);
ix = indsT + (indsC-1) * nT;

% extract only spatial indices within the col index
clips = reshape(chrawdat(ix), numel(dt), numel(dc), numel(row));
electwaveforms = squeeze(clips(:,:));

% normelecwf= zeros(size(electwaveforms));
% for ii = 1:size(electwaveforms,2)
%     normelecwf(:,ii) = electwaveforms(:,ii) / max(electwaveforms(:,ii));
% end
%
% CC = CC + (c * c')/1e3;
% [U Sv V] = svdecon(CC);
% wPCA = U(:, 1:3);
% wPCA(:,1) = - wPCA(:,1) * sign(wPCA(8,1));
%

% replace pulses with median of the signal (generally zero)
vtoreplace = median(chrawdat,'all');
chrawdat(squeeze(ix)) = vtoreplace;

%mu = ones(size(row));
%for n=1:60; plot(1:size(chrawdat,1),chrawdat(:,n),row(col==n),mu(col==n),'ro'); title(n); pause; end

chdat = int16(chrawdat');
% get electric pulses and thier time per channel
electpulsedat   = arrayfun(@(x)(electwaveforms(:,col==x)'),1:nChan,'un',0)';
% get time of each  electric event
electpulsetimes = chtvec(squeeze(indsT));
pulsestartstop  = arrayfun(@(x)(electpulsetimes(:,col==x)'),1:nChan,'un',0)';

fprintf('Took %2.2f seconds...\n',toc);



end


function [chdat, pulsestartstop, electpulsedat] = remelecpulses_old(eventstream, dat, chanMap)


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

madTh           =   30;
sumsig          =   abs(sum(chdat,1)); % take absolute sum of all signals
Thres           =   madTh * mad(diff(sumsig));
fonsets         =   find( diff(diff(sumsig) >  Thres) <0 ) + 1;

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


% p12 = (chtvec > pulsestart(ii)-(1.2*tbefpulse(1)) & chtvec < pulsestart(ii+1));
%
% close all
% plot(dat.ChannelDataTimeStamps(p12),chdat(12,p12),pulsestart(ii),300,'ro')
% hold on
% plot(dat.ChannelDataTimeStamps(tbefpulse),chdat(12,tbefpulse),'k','LineWidth',3)
% plot(dat.ChannelDataTimeStamps(taftpulse),chdat(12,taftpulse),'m','LineWidth',3)
%n = 7; plot(chtvec,chdat(n,:),'-',stimstrtpts,0,'rx',stimendpts,0,'ko');

end




function S1 = pulseminpoint(S1, sig, varargin)
% takes an extra argument which specifies which dimension to filter on
% extra argument can be a vector with all dimensions that need to be
% smoothed, in which case sig can also be a vector of different smoothing
% constants

idims = 2;
if ~isempty(varargin)
    idims = varargin{1};
end
if numel(idims)>1 && numel(sig)>1
    sigall = sig;
else
    sigall = repmat(sig, numel(idims), 1);
end

for i = 1:length(idims)
    sig = sigall(i);
    
    idim = idims(i);
    Nd = ndims(S1);
    
    S1 = permute(S1, [idim 1:idim-1 idim+1:Nd]);
    
    dsnew = size(S1);
    
    S1 = reshape(S1, size(S1,1), []);
    dsnew2 = size(S1);
    
    S1 = cat(1, Inf*ones([sig, dsnew2(2)]),S1, Inf*ones([sig, dsnew2(2)]));
    Smax = S1(1:dsnew2(1), :);
    for j = 1:2*sig
        Smax = min(Smax, S1(j + (1:dsnew2(1)), :));
    end
    
    S1 = reshape(Smax, dsnew);
    
    S1 = permute(S1, [2:idim 1 idim+1:Nd]);
end
end



