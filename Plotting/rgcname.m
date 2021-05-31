

function [rgctitlename, pngfilename] = rgcname(stimname, chinfo, experimentdate, iter, varargin)
%
%%% rgcname %%%
%
% This function generate string name for each cell based on the format that
% is used for naming analyzed data of each cells.
%
%
% ===============================Inputs====================================
%
%   stimname : name of the stimulus pluse onther necessary texts.
%   chinfo : channel info from the excel file generate from spike sorting.
%   datapath : folder path of the data.
%
%================================Output====================================
%
%   rgcname : output RGC name.
%   rgcpathname : folder path for the input RGC to be saved.
%
% written by Mohammad, 05.12.2016.

if nargin < 3 || isempty(experimentdate)
    experimentdate = date;
end

if nargin > 4
    breakto2line = varargin{1};
else
    breakto2line = false;
end

if isstruct(chinfo)
    if breakto2line
        rgctitlename ={[stimname,' Analysis for Channel ',num2str(chinfo.ch),', Cluster ',num2str(chinfo.clus),...
            ', KS ID ',num2str(chinfo.id)];['quality ',num2str(chinfo.quality),', nspk ',num2str(chinfo.n_spikes),...
            ', ',strrep(chinfo.comment{1},'_',' '),', for Experiment on ',experimentdate]};
        pngfilename = [rgctitlename{1},', quality ',num2str(chinfo.quality), ', expdate ',experimentdate];
    else
        rgctitlename =[stimname,' Analysis for Channel ',num2str(chinfo.ch),', Cluster ',num2str(chinfo.clus),...
            ', KS ID ',num2str(chinfo.id),', quality ',num2str(chinfo.quality),', nspk ',num2str(chinfo.n_spikes),...
            ', ',strrep(chinfo.comment{1},'_',' '),', for Experiment on ',experimentdate];
        pngfilename = [extractBefore(rgctitlename,', nspk '),', expdate ',experimentdate];
    end
    
else
    
    if breakto2line
        rgctitlename ={[stimname,' Analysis for Channel ',num2str(chinfo(1,1)),', Cluster ',num2str(chinfo(1,2)),','];...
            [' for Experiment on ',experimentdate]};
        pngfilename = [rgctitlename{1},rgctitlename{2}];
    else
        rgctitlename =[stimname,' Analysis for Channel ',num2str(chinfo(1,1)),', Cluster ',num2str(chinfo(1,2)),...
            ', for Experiment on ',experimentdate];
        pngfilename = rgctitlename;
    end
    
end

rgctitlename = [num2str(iter,'%02g-'), rgctitlename];
pngfilename = [num2str(iter,'%02g-'), pngfilename];

end