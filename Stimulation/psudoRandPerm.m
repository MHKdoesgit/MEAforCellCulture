
function varargout = psudoRandPerm(seed,inputVec,varargin)
%
%%% psudoRandPerm %%%
%
%
% This function generate psudorandom permution similar to randperm in MATLAB
% but works with psudorandom number generator ran1. It also gives back the
% most recent seed value to continue the permuation in case of repeated trials.
% note that the direction of permutation is along x-axix or for columns of
% MATALB not for the rows.
%
%
% ===============================Inputs====================================
%
%   seed : seed value for random number generation.
%   inputVec : input vector used for permutation.
%
%================================Output====================================
%
%   testOrder : vector of permuted indices for the inputVec.
%   newSeed : the recent seed that used in the ran1 function.
%   outputVec : the permuted input vector along x-axis
%
% Note that the permution algorithem is based on Fisher-Yater shuffle
% algorithem identical to what is used in the stimulus program.
% for more info check :
% https://en.wikipedia.org/wiki/Fisher%E2%80%93Yates_shuffle
%
% written by Mohammad, 01.02.2016

newSeed = seed;

testOrder = zeros(1,size(inputVec,2));
testOrder(1) = 1;

for i = 2:length(inputVec)-1
    
    [randVal,newSeed] = ran1(newSeed);
    
    j = ceil(i*randVal);    % based on Fischer-Yates algorithem
    testOrder(i) = testOrder(j);
    testOrder(j) = i;
    
end

testOrder= [testOrder(end),testOrder(1:end-1)]+1;   % to match MATLAB indexing
varargout{1} = testOrder;
varargout{2} = newSeed;
for j = 1:size(inputVec,1)
    varargout{3}(j,:) = inputVec(j,testOrder);
end;
end
