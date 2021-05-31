

function [resplate, resppeak, resplatefit, resppeakfit] = psthPeakLatency(psth, activreg, varargin)
%
%%% psthPeakLatency %%%
%
%
% This function finds the peak of the input psth and measures the
% time-to-peak value of the input psth.
%
%================================Inputs====================================
%
%   psth : input psth.
%   activreg : time-range in which the peak should be measured.
%   thresh : the threshold to ignore the very small peak.
%   stimdur : duration of the stimulus (by default 500 ms).
%
%================================Output====================================
%
%   resplate : peak latency at the time resolution of the input psth.
%   resppeak : psth peak.
%   resplatefit : psth latency of the fit at the location of the peak.
%   resppeakfit : peak values of the fit at the resppeak location.
%
% written by Mohammad, 06.08.2017.

if nargin > 2, thresh = varargin{1}; else, thresh = 0; end
if nargin > 3, stimdur = varargin{2}; else, stimdur = max(size(psth)); end

if (length(activreg)==1 ) && (activreg > 0)
    activreg = [25, 250];
elseif (length(activreg)==1 ) && (~activreg)
    activreg = [0, stimdur];
end
warning('off','all');
%if all(activreg) < 10,    activreg = activreg *100;     end;
%if stimdur < 10, stimdur = stimdur*100; end;

% tvec = linspace(0, stimdur, length(1:0.01:size(psth,1)));
tvec = linspace(0, stimdur, size(psth,1));
activreghr = find(tvec >= activreg(1) & tvec <= activreg(end));
[resplate, resppeak, resplatefit, resppeakfit]  = deal(nan(1, size(psth,2)));

%smfun = @(x)(interp1(1:size(psth,1), x, 1:0.01:size(psth,1),'pchip'));

for ii = 1: size(psth,2)
    
    [peakval, peakidx] = findpeaks((psth(:,ii)),'npeaks',5,'sortstr','descend','minpeakdistance',1);
    if length(peakidx)< 5
        peakidx = [peakidx;repmat(max(activreghr),5-length(peakidx),1)]; %#ok
        peakval = [peakval;zeros(5-length(peakval),1)];     %#ok
    end
    
    if ismember(peakidx(1),activreghr)
        finalpeakidx = peakidx(1);
        peakmagidx = 1;
        
    elseif ismember(peakidx(2),activreghr)
        finalpeakidx = peakidx(2);
        peakmagidx = 2;
        
    elseif ismember(peakidx(3),activreghr)
        finalpeakidx = peakidx(3);
        peakmagidx = 3;
        
    elseif ismember(peakidx(4),activreghr)
        finalpeakidx = peakidx(4);
        peakmagidx = 4;
        
        
    elseif ismember(peakidx(5),activreghr)
        finalpeakidx = peakidx(5);
        peakmagidx = 5;
        
    else
        [~,actregidx] = findClosestValue(peakidx, max(activreghr));
        finalpeakidx = peakidx(actregidx);
        peakmagidx = actregidx;
    end
    
    % now double check if there was a bit smaller peak just before
    if (peakval(1) - peakval(2) <= 0.4*peakval(1)) && ismember(peakidx(2),activreghr) && ...
            (peakidx(2) < peakidx(1))
        finalpeakidx = peakidx(2);
        peakmagidx = 1;
    end
    
    if (peakval(1) - peakval(3) <= 0.4*peakval(1)) && ismember(peakidx(3),activreghr) ...
            && (peakidx(2) > peakidx(3)) && (peakidx(3) < peakidx(1))
        finalpeakidx = peakidx(3);
        peakmagidx = 1;
    end
    
    
    if (peakval(1) - peakval(4) <= 0.4*peakval(1)) && ismember(peakidx(4),activreghr) ...
            && (peakidx(2) > peakidx(4)) && (peakidx(3) > peakidx(4)) && (peakidx(4) < peakidx(1))
        finalpeakidx = peakidx(4);
        peakmagidx = 1;
    end
    
    if (peakval(1) - peakval(5) <= 0.4*peakval(1)) && ismember(peakidx(5),activreghr) ...
            && (peakidx(2) > peakidx(5)) && (peakidx(3) > peakidx(5)) && (peakidx(4) > peakidx(5)) ...
            && (peakidx(5) < peakidx(1))
        finalpeakidx = peakidx(5);
        peakmagidx = 1;
    end
    
    if finalpeakidx == size(psth,1)
        finalpeakidx = finalpeakidx-1;
    end
    
    resplate(ii) = tvec(finalpeakidx);
    resppeak(ii) = peakval(peakmagidx);
    
    % small polynomial fit at the peak to get more precise latency
    tvechr = linspace(tvec(finalpeakidx-1), tvec(finalpeakidx+1),1e3);
    pfit = polyval(polyfit(tvec(finalpeakidx-1:finalpeakidx+1),psth(finalpeakidx-1:finalpeakidx+1,ii)',2), tvechr);
    [resppeakfit(ii),pfitidx] = max(pfit);
    resplatefit(ii) = tvechr(pfitidx);
    
    resplate(ii) = correctNonresppeak(resppeak(ii), resplate(ii), thresh, max(activreg));
    resplatefit(ii) = correctNonresppeak(resppeak(ii), resplatefit(ii), thresh, max(activreg));
    
    clearvars peakval peakidx tvechr pfit pfitidx  finalpeakidx peakmagidx;
end
warning('on','all');

end

%--------------------------------------------------------------------------------------------------%
%----------                               sub-functions                                  ----------%
%--------------------------------------------------------------------------------------------------%

function lateout = correctNonresppeak(peak, late,threshold, maxlate)

if peak < threshold
    lateout = maxlate;
else
    lateout = late;
end

end
