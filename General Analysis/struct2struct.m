

function s_output = struct2struct (s_input, s_tojoin, varargin)
%
%%% struct2struct %%%
%
%
% This function copies sub-fields of a structure to new structure with same
% order of sub-fields. It reduces one-dimention in copying structures.
%
% ===============================Inputs====================================
%
%   s_input : input structure.
%   s_tojoin : second structure to join to the first input.
%   flipflag : flag to put subfield of second input on top of first input.
%
%================================Output====================================
%
%   s_output : out-put structure.
%
% written by Mohammad, 13.11.2015.
% added option for array structures and flag for flipping the order on
% 07.02.2017 by Mohammad.

% check inputs
if nargin > 2
    flipflag = logical(varargin{1});
else
    flipflag = false;
end

% to put the second input fields on top of first input
if flipflag
    s_output = s_tojoin;
    s_tojoin = s_input;
else
    s_output = s_input;
end

s_fields = fieldnames(s_tojoin);
% for array structures
for ii = 1:length(s_input)
    for jj = 1:numel(s_fields)
        s_output(ii).(s_fields{jj}) = s_tojoin(ii).(s_fields{jj,:});
    end
end

end