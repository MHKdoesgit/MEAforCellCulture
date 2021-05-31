

function varargout = Retinalabeler (datapath,varargin)
%
%%% retinaLabeler %%%
%
%
% This function create long name of the experiment folder based on the
% short name format and I use generally.
%
%
% ===============================Inputs====================================
%
%   datapath : folder path of the experiment.
%
%================================Output====================================
%
%   retinaLongName : the long name format of the experiment folder.
%
% written by Mohammad, 28.04.2015
% update to completely new version with retina side option on 13.05.2016

retNameLoc = regexp(datapath,{'fr_fp','fr_sp','fr_tp','sr_fp','sr_sp','sr_tp'});

if sum(not(cellfun(@isempty,retNameLoc))) == 0
    error('the input name does not contain any of the valid retina names,check and come back bro!');
elseif sum(not(cellfun(@isempty,retNameLoc))) == 1
    slashpos = max(strfind(datapath,'\'));
    if slashpos > cell2mat(retNameLoc)
        retName = datapath(cell2mat(retNameLoc):slashpos);
    else
        retName = datapath(cell2mat(retNameLoc):end);
    end
end

switch lower(retName)
    case {'fr_fp','fr_fp_lv','fr_fp_ld','fr_fp_rv', 'fr_fp_rd'}
        retinaName = 'first retina, first piece';
    case {'fr_sp','fr_sp_lv','fr_sp_ld','fr_sp_rv', 'fr_sp_rd'}
        retinaName = 'first retina, second piece';
    case {'fr_tp','fr_tp_lv','fr_tp_ld','fr_tp_rv', 'fr_tp_rd'}
        retinaName = 'first retina, third piece';
    case {'sr_fp','sr_fp_lv','sr_fp_ld','sr_fp_rv', 'sr_fp_rd'}
        retinaName = 'second retina, first piece';
    case {'sr_sp','sr_sp_lv','sr_sp_ld','sr_sp_rv', 'sr_sp_rd'}
        retinaName = 'second retina, second piece';
    case {'sr_tp','sr_tp_lv','sr_tp_ld','sr_tp_rv', 'sr_tp_rd'}
        retinaName = 'second retina, third piece'; 
end

retside = retName(end-1:end);
switch lower(retside)
    case 'lv'
        retinaSide = 'left eye, ventral retina';
    case 'ld'
        retinaSide = 'left eye, dorsal retina';
    case 'rv'
        retinaSide = 'right eye, ventral retina';
    case 'rd'
        retinaSide = 'right eye, dorsal retina';
    otherwise
        retinaSide = [];
end
        
retinaLongName = [retinaName,', ',retinaSide];        
varargout{1} = retinaLongName;

end
