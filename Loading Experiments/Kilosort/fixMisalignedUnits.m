
function [spktimes, avgtemps, misalUnitIds] = fixMisalignedUnits(spktimes, spkclusts, spktemps, temps, comments)
% Returns fixed spike times and averaged templates for peak channel extraction
% Inputs:
% spktimes, spkclusts, spktemps correspond to KS output (all same size)
% spkclusts should run from 1 to Nunits
% temps are the templates from KS
% comments is a cell array Nunits x 1 with the text comments
% Outputs:
%       spktimes: for all units, misaligned units are fixed
%       avgtemps: average template for all units, usefull for finding peak channel
%==========================================================================
[Ntemps, Nt, Nyx] = size(temps);
Nunits = max(spkclusts);
temps (isnan(temps(:))) = 0; % why the values come as NaN from KS
%==========================================================================
% get average templates
templateMat = accumarray([spkclusts spktemps+1], 1, [Nunits Ntemps]);

avgtemps = (templateMat * temps(:,:))./sum(templateMat, 2);
avgtemps = reshape(avgtemps, [Nunits Nt Nyx]);
%==========================================================================
% fix misalignment in marked units (do it for all in future)
fprintf('Fixing misalinged spikes...'); tic; msg = [];
if iscell(comments)
    if all(contains(unique(comments),{'good','noise'})) % get good cells
        misalUnitIds = find(contains(comments, 'good'));
    else
        misalUnitIds = find(contains(comments, 'misaligned'));
    end
elseif ~iscell(comments) && strcmpi(comments,'all')
    misalUnitIds =  unique(spkclusts);
end
maxLag = ceil(size(temps, 2) / 2);

for ii = 1:numel(misalUnitIds)
    
    unitid = misalUnitIds(ii);
    spkids = spkclusts == unitid;
   
    unittimes = spktimes(spkids);
    unittemps = spktemps(spkids)+ 1;
    
    tempids = unique(unittemps);
    ctemps = temps(tempids,:,:);
    
    [~, bestid] = min(min(ctemps(:,:), [], 2));
    finaltemp = zeros(Nt, Nyx, 'single');
    temp1 = squeeze(ctemps(bestid,:,:));
    for it = 1:numel(tempids)
        temp2 = squeeze(ctemps(it,:,:));
        [~, imax]= maxSlidingCorr(temp1', temp2', maxLag);
        temptimeids = unittemps == tempids(it);
        unittimes(temptimeids) = unittimes(temptimeids) - imax;
        if imax>=0
            finaltemp(imax+1:Nt, :) = finaltemp(imax+1:Nt, :) + nnz(temptimeids) * temp2(1:(Nt-imax), :);
        else
            finaltemp(1:(Nt + imax), :) = finaltemp(1:(Nt + imax), :) + nnz(temptimeids) * temp2(-imax+1:Nt, :);
        end
    end
    finaltemp = finaltemp/numel(unittimes);
    avgtemps(unitid, :, :) = finaltemp;
    
    spktimes(spkids) = unittimes;
    %--------------------------------------------------------------------------
    fprintf(repmat('\b', 1, numel(msg)));
    msg = sprintf('Units %d/%d. Time elapsed %2.2f s...\n', ii, numel(misalUnitIds),toc);
    fprintf(msg);
end

end

function [mcorr, imax]= maxSlidingCorr(temp1, temp2, maxLag)

[Nyx, Nt]= size(temp1);

lag = -maxLag:maxLag;
XCmat = zeros(2, 2, length(lag), 'single');
spike_templates = reshape(cat(3,temp1, temp2), Nyx, Nt, 2);

for ilag = 0:maxLag
    %shift matrices

    Y1 = reshape(spike_templates(:, 1+ilag:end, :), Nyx * (Nt - ilag), 2); 
    Y2 = reshape(spike_templates(:, 1:end-ilag, :), Nyx * (Nt - ilag), 2);
    resmat = corr(Y1,Y2); %core calculation
    
    XCmat(:,:, maxLag+1+ilag) = resmat; 
    XCmat(:,:, maxLag+1-ilag) = resmat';
end

[mcorr, im] = max(XCmat(1,2,:));
imax = lag(im);

end

