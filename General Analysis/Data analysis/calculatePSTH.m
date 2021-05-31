

function [ allRates, myTimes] = calculatePSTH( trialRaster, timeLims, varargin)
%CALCULATEPSTH Summary of this function goes here
%   Input:
%       trialRaster: cell array, cellN x trialN
%       dimInd: index of which dimension different trial are (default is last)
%   Output:
%       allRates: 
%Impovements:
% 2) give SD in this function as well with varargout (don't calculate if not needed)


%timestep as third input
if (nargin > 2 && isnumeric(varargin{1})); timeStep = varargin{1};  else, timeStep= 0.01; end
%position of averaging dimension (trials)
if nargin > 3 && isinteger(varargin{2}); dimInd = varargin{2}; else, dimInd = ndims(trialRaster); end
%==========================================================================
%Create edges
n=floor(diff(timeLims)/timeStep);
myEdges=timeLims(1):timeStep:timeLims(1)+n*timeStep;
if n*timeStep <timeLims(2)
    myEdges=[myEdges timeLims(2)];
end
myTimes=myEdges(1:end-1)+diff(myEdges)/2;
%==========================================================================

targetSize=size(trialRaster); targetSize(dimInd)=[];
allRates=zeros(numel(myTimes), prod(targetSize));

trialN=size(trialRaster,dimInd);
trialRaster=permute(trialRaster,[dimInd 1:dimInd-1]);

trialRaster=reshape(trialRaster, [trialN prod(targetSize)]);

C=cellfun(@(x) sum(isnan(x)), trialRaster);
trialNall=trialN-sum(C,1);

for ii=1:size(trialRaster,2)
     myMat=CelltoMatUE(trialRaster(:,ii));
     myRates=histc(myMat(:), myEdges)'; myRates(end)=[]; %counting is happening here
     allRates(:,ii)=myRates(:);
end
%use bsxfun to divide all of them at the same time
allRates=bsxfun(@times,allRates,1./trialNall);
divelement=1./diff(myEdges(:));
allRates=bsxfun(@times,allRates,divelement);
allRates=reshape(allRates,[numel(myTimes) targetSize]);
allRates=permute(allRates, [2:ndims(allRates) 1]);
%targetSize(targetSize~=numel(myTimes))=1; divelement=zeros(targetSize);
%divelement(:)=1./((trialN)*diff(myEdges(:))');
%allRates=bsxfun(@times,allRates,divelement);
end

