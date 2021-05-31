

function s_out = structarray2arraystruct(s_in, varargin)
%
%%% structarray2arraystruct %%%
%
%
% This function removes one layer of multi-layer array structure and turn
% it into a array structure with one column. It is equivalant to (:) for
% the vectors, but here for structures.
%
% ===============================Inputs====================================
%
%   s_input : multi layer input structure.
%
%================================Output====================================
%
%   s_output : out-put structure.
%
% written by Mohammad, 07.11.2015.

s_in_fields = fieldnames(s_in);
s_in_maxsize = max(structfun(@numel,s_in));

temp_cell = cell(s_in_maxsize,numel(s_in_fields));

for ii = 1: max(size(fieldnames(s_in)))
    
    for jj = 1:length(s_in.(s_in_fields{ii}))
        temp_cell{jj,ii} = s_in.(s_in_fields{ii})(jj);
    end
    
end
s_out = temp_cell(:);
s_out = cell2mat(s_out(not(cellfun('isempty',temp_cell(:)))));

end