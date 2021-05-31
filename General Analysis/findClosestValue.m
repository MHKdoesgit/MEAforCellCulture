
function [varargout] = findClosestValue (x,values, varargin)
%
%%% findClosestValue %%%
%
%
% This function find the closest value in the vector x to the input value.
% if a second vector is also used, it finds the value of first vetor inside
% second vector.
%
%================================Inputs====================================
%
%    x : input vector.
%    value : target value.
%    varargin{1} : optional, second vector (y).
%
%================================Output====================================
%
%   varargout{1} : closest value of x or y to the target value.
%   varargout{2} : index of the target value.
%
% written by Mohammad, 24.05.2015
% added option for vector values 05.06.2015

Idx = zeros(1,length(values));
out = zeros(1,length(values));

for ii = 1:length(values)
    [~,Idx(ii)] = min(abs(x - values(ii)));
    
    if (nargin < 3)
        out(ii) = x(Idx(ii));
    else
        y = varargin{1};
        out(ii)= y(Idx(ii));
    end
end

varargout{1} = out;
varargout{2} = Idx;