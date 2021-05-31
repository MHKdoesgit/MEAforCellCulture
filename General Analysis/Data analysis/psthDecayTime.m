
function [decayidx, decayfit, decaytvec] = psthDecayTime(psth, peaklatency, varargin)
%
%%% psthDecayTime %%%
%
%
% This function calculates the decay time of the input PSTH.
%
% ===============================Inputs====================================
%
%   psth : input psth.
%   peaklatency : time to peak, helps for peak detection. to get it use
%                 psthpeaklatency function.
%   stimdur : duration of the stimulus in ms (defualt is 500).
%
%================================Output====================================
%
%   decayidx : decay time for the input psth.
%   decayfit : exponential fit that was used to calculate the decay time.
%   decaytvec : time vector for the exonential fit.
%
% written by Mohammad, 08.11.2018.

if nargin > 2, stimdur = varargin{1}; else, stimdur = 500; end

tvec = linspace(0,stimdur,size(psth,1));
[decayidx, sparesp] = deal(zeros(1, size(psth,2)));
[decayfit, decaytvec] = deal(cell(1, size(psth,2)));

for ii = 1: size(psth,2)
    
    [~,peakidx] =  findClosestValue(tvec, peaklatency(ii));
    
    [param, decayfit{ii}] = exponential_fit(tvec(peakidx:end), psth(peakidx:end,ii)');
    decaytvec{ii} = tvec(peakidx:end);
    
    %fitvals = exp(polyval(p,tvec(peakidx:end)));
    decayidx(ii) = (-1/ param(2))*1e3; % in milisecond
    sparesp(ii) = param(3);
    
end

if size(decayfit,1) == size(decayfit,2)
    decayfit = decayfit{1};
    decaytvec = decaytvec{1};
end

end

%--------------------------------------------------------------------------------------------------%

function [param, ypred, f] = exponential_fit(x,y)

x = (x- x(1))/1e3;  % first bring it to zero and then in second
f = @(B,x)  B(1).*exp(B(2).*x) + B(3);     % B(1) = a,  B(2) = b,  B(3) = c

B0 = [y(1) -100 min(y)];  % % B0(1) = peak firing rate,
% B0(2) = estimate decay time in second based on alhpa cells paper from Meister (2017) -1/0.1
% B0(3) = last point or spontanous firing rate
warning('off','all');
try
    [BETA,RESID,J,COVB,MSE] = nlinfit(x,y,f,B0);    %#ok
catch 
    %disp(ME.message);
    BETA = nan(1,3);
end
warning('on','all');
param = BETA;
ypred = f(BETA,x);

end