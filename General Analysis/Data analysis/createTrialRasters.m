function [ allRasters ] = createTrialRasters(spikeTimes, onTimes, offTimes, varargin)
%CREATETRIALRASTERS Summary of this function goes here
%   Detailed explanation goes here
%make more flexible if on and offtimes have NaNs - added, shouldn't hurt

assert(isequal(size(onTimes), size(offTimes)), 'onTimes and offTimes must have the same size!');

if nargin<4; shift=0; else, shift=varargin{1}; end

cellN=numel(spikeTimes);
extraDimensions=size(onTimes); extraDimensions(extraDimensions==1)=[];%fix when using 1 cell
trialN=prod(extraDimensions);
allRasters=cell(trialN,cellN);

for cellId=1:cellN
    spk=spikeTimes{cellId};
    for ii=1:trialN
        ons=onTimes(ii); offs=offTimes(ii);
        if or(isnan(ons),isnan(offs))
            val=NaN;
        else 
            val=spk(spk>ons & spk<offs)-ons-shift;
        end
        allRasters{ii,cellId}=val;
    end
end
allRasters=reshape(allRasters', [cellN extraDimensions]);

%
% for ii=1:numel(allRasters)
%      [outputs{:}] =ind2sub(finalDimensions, ii);
%      spk=spikeTimes{outputs{1}};
%      ons=onTimes(outputs{2:end}); offs=offTimes(outputs{2:end});
%      if or(isnan(ons),isnan(offs)); val=NaN;
%      else val=spk(spk>ons & spk<offs)-ons-shift; end;
%      allRasters{outputs{:}}=val;
% end

%old code - might be a bit faster, not general across dimensions
% cellN=numel(spikeTimes);
% trialN=numel(onTimes);
%
% allRasters=cell(cellN, trialN);
%
% for cellId=1:cellN
%     spk=spikeTimes{cellId};
%     for ii=1:trialN
%         trialTimes=spk(spk>onTimes(ii) & spk<offTimes(ii));
%         allRasters{cellId, ii}=trialTimes(:)'-onTimes(ii);
%     end
% end


end

