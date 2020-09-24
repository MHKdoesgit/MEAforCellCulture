
function normalized = normalizetoRange(inputArray, minval, maxval)
%
%%% normalizetoRange %%%
%
%
% This function normalize the input array to the range defined by X and Y
%
% ===============================Inputs====================================
%
%   inputArray : input array.
%   minval : minimum range for normalization.
%   maxval : maximum range for normalization.
%
%================================Output====================================
%
%   normalized : normalized input array in range of X and Y.
%
% written by Mohammad, 21.10.2015, based on Max function at 
% http://stackoverflow.com/questions/10364575/normalization-in-variable-range-x-y-in-matlab

     % Normalize to [0, 1]:
     size_array = size(inputArray);
     inputArray = inputArray(:);
     m = min(inputArray);
     range = max(inputArray) - m;
     inputArray = (inputArray - m) / range;

     % Then scale to [x,y]:
     range2 = maxval - minval;
     normalized = reshape((inputArray*range2) + minval,size_array);
     
end