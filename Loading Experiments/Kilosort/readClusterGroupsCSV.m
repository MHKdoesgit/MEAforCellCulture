

function [cgs, cids] = readClusterGroupsCSV(filename)
%function [cids, cgs] = readClusterGroupsCSV(filename)
% cids is length nClusters, the cluster ID numbers
% cgs is length nClusters, the "cluster group":
% - 0 = noise
% - 1 = mua
% - 2 = good
% - 3 = unsorted

% update to completely new version by Mohammad on 19.02.2020.

headerseq = getTSVheader(filename);
fid = fopen(filename);
% C = textscan(fid, '%s', 'Delimiter', '\t','EmptyValue', Inf);
C = textscan(fid, repmat('%s',1,numel(headerseq)), 'Delimiter', '\t');
fclose(fid);

for ii = 1: numel(headerseq) % fix this later
    phyclusinfo = cellfun(@str2num, C{ii}(2:end), 'uni', false);
    if all(cellfun('isempty',phyclusinfo))
        phyclusinfo = C{ii}(2:end);
        emptyIndex = cellfun('isempty', phyclusinfo);     % Find indices of empty cells
        switch lower(C{ii}{1})
            case 'kslabel'
                phyclusinfo(emptyIndex) = {'phyclusters'};
            case 'group'
                if any(emptyIndex)
                    warning('Aint nobody sorted the channels below, group label set to: unsorted');
                    disp(repmat('-',1,50)); disp('          id         phy-ch');     disp(repmat('-',1,50));
                    disp([cgs.id(emptyIndex),cgs.ch(emptyIndex)]);
                    phyclusinfo(emptyIndex) = {'unsorted'};
                end
        end
        
    else
        emptyIndex = cellfun('isempty', phyclusinfo);     % Find indices of empty cells
        switch lower(C{ii}{1})
            case {'quality','qual','q'}
                phyclusinfo(emptyIndex) = {max(cell2mat(phyclusinfo))+1};
            case {'ch', 'channel'}
                phyclusinfo = cellfun(@(x)x+1,phyclusinfo,'un',0); % add 1 to match indexing from MATLAB
            otherwise
                phyclusinfo(emptyIndex) = {NaN};
        end
        phyclusinfo = single(cell2mat(phyclusinfo));
    end
    cgs.(C{ii}{1}) =  phyclusinfo;
end

if numel(headerseq)<=2
    cids = round(single((cgs.cluster_id)));
    propval = cgs.(C{end}{1});
    
    if iscell(propval)
        cgs = [num2cell(cids),propval];
    else
        cgs = [cids, single(propval)];
    end
else
    cids = cgs.id;
end

% this is the old code, from the KS/phy peoples, it's good for shit!
% 
% cids = cellfun(@str2num, C{1}(2:end), 'uni', false);
% ise = cellfun(@isempty, cids);
% cids = [cids{~ise}];
% 
% isUns = cellfun(@(x)strcmp(x,'unsorted'),C{2}(2:end));
% isMUA = cellfun(@(x)strcmp(x,'mua'),C{2}(2:end));
% isGood = cellfun(@(x)strcmp(x,'good'),C{2}(2:end));
% 
% cgs = zeros(size(cids));
% 
% cgs(isMUA) = 1;
% cgs(isGood) = 2;
% cgs(isUns) = 3;

end

%--------------------------------------------------------------------------------------------------%
%----------                               sub-functions                                  ----------%
%--------------------------------------------------------------------------------------------------%

function headerseq = getTSVheader(filename)

fid = fopen(filename);
% C = textscan(fid, '%s%s');
D = textscan(fid, '%s','Delimiter', '\t');
fclose(fid);
n = cellfun(@(x)(~isnan(str2double(x))),D{1});
headerseq = D{1}(1:find(n,1,'first')-1);

end
