


function st_out = struct2subfield(st_in,newfield,field2sub,varargin)
%
%%% struct2subfield %%%
%
%
% This function get the names of some of the fields of the input structure
% and put them inside one new structure subfield with user defined name.
% It is mostly use to organize structures.
%
%================================Inputs====================================
%
%   st_in : input structure.
%   newfield : name of the new field.
%   field2sub : names of fields to move to new field.
%
%================================Output====================================
%
%   st_out : output structure with sub-field structure.
%
% written by Mohammad, 10.02.2017.

if nargin > 3, 
    rmidx = varargin{1}:numel(field2sub); 
else
    rmidx = 1:numel(field2sub);
end;

if (nargin > 3 && length(varargin{1})>1)
    rmidx = varargin{1};
end;

if (nargin > 4 )
    rmidx = varargin{1}:varargin{2};
end;

sttemp = st_in;
st_out = cell(size(st_in));

for i= 1:size(st_in,1)
    for j = 1:numel(field2sub)
        sttemp(i).(newfield).(field2sub{j}) = st_in(i).(field2sub{j});
    end
    st_out{i} = rmfield(sttemp(i),field2sub(rmidx));
end;
st_out = cell2mat(st_out);


end