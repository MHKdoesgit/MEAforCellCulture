

function ccdata = calcCCG(spiketimes, timebinsms, lagnumbins, varargin)
%
%%% calcCCG %%%
%
%
% This function calculates cross-correlogram from all the spikes included
% in the spike-time. In this function, instead of binning the spikes the
% region around every spikes is cheked for the auto and cross correlogram.
% This function is comes from very similar function offered as part of
% Kilosort 2 software.
%
%================================Inputs====================================
%
%   spiketimes : spike times.
%   ftimes : frametimes, used to calculate experiment duration (not
%            necessary).
%   timebinsms : time bins for the correlogram in ms (e.g. 0.5 or 0.25).
%   lagnumbins : number of bins for the lag around the peak (e.g. 50 or 100).
%   batchlen : for datasets with many cells, this can work as batch number.
%
%================================Output====================================
%
%   ccdata : correlogram data, includes cross-correlogram,
%            auto-correlogram, timebins etc.
%
% written by Mohammad, 05.05.2020.

if nargin > 4, batchlen = varargin{1}; batchflag = true; else, batchflag = false; end

% this function merges clusters based on template correlation
% however, a merge is veto-ed if refractory period violations are introduced

%ops = rez.ops;
%dt = 0.5/1000;
lagamount = floor(lagnumbins/1e3/timebinsms);

% Xsim = rez.simScore; % this is the pairwise similarity score
% Nk = size(Xsim,1);
% Xsim = Xsim - diag(diag(Xsim)); % remove the diagonal of ones

Nk = numel(spiketimes);
if batchflag, Nkiter = batchlen; else, Nkiter = Nk;   end
% if nargin < 5
%     Nkiter = Nk;    
% else
%     Nkiter = batchlen;%varargin{1};
% end
if numel(Nkiter)==1, Nkiter = [1,Nkiter]; end
% sort by firing rate first
ccgmetadata.numspk = cellfun('length',spiketimes(Nkiter(1):Nkiter(2)))';
% since there is no ftimes here, we measure the experiment duration as
% difference between the first and last spike
ccgmetadata.expdur = max(cell2mat(spiketimes))- min(cell2mat(spiketimes));
ccgmetadata.frate = (ccgmetadata.numspk ./ ccgmetadata.expdur)';
% get location info
% ccgmetadata.meanspkloc = CelltoMatUE(cellfun(@nanmean,loc(Nkiter(1):Nkiter(2)),'un',0));
% ccgmetadata.stdspkloc = CelltoMatUE(cellfun(@nanstd,loc(Nkiter(1):Nkiter(2)),'un',0));
ccgmetadata.tbin = timebinsms;
%ccgmetadata.lagbins = lagamount;
ccgmetadata.lag = timebinsms*(-lagamount:lagamount);
% nspk = zeros(Nk, 1);
% for j = 1:Nk
%     nspk(j) = sum(rez.st3(:,2)==j); % determine total number of spikes in each neuron
% end
%[~, isort] = sort(nspk,'descend'); % we traverse the set of neurons in ascending order of firing rates
% fprintf('initialized spike counts\n')

% pre-allocating memory
cc = zeros(Nk,length(Nkiter(1):Nkiter(2)),2*lagamount+1,'single');
[qc, rc] = deal(zeros(Nk,length(Nkiter(1):Nkiter(2)),'single'));

% tic;
iter = 1;
for j = Nkiter(1):Nkiter(2)
%     s1 = spk{isort(j)}; %rez.st3(rez.st3(:,2)==isort(j), 1)/ops.fs; % find all spikes from this cluster
%     if numel(s1)~=nspk(isort(j))
%         fprintf('lost track of spike counts') %this is a check for myself to make sure new cluster are combined correctly into bigger clusters
%     end
    s1 = spiketimes{j};
    %ienu = 1:Nk;

    K = zeros(Nk,2*lagamount+1,'single');
    [Q,R]= deal(zeros(Nk,1,'single')); 
    
    parfor k = 1:Nk
        s2 = spiketimes{k}; % find the spikes of the pair
        [K(k,:), Qi, Q00, Q01, rir] = ccgfast(s1, s2, lagamount, timebinsms);
        Q(k) = min(Qi/(max(Q00, Q01))); % normalize the central cross-correlogram bin by its shoulders OR by its mean firing rate
        R(k) = min(rir); % R is the estimated probability that any of the center bins are refractory, and kicks in when there are very few spikes
        
        
    end
    K = K ./ max(K(j,:));

    cc(:,iter,:) = K;
    qc(:,iter) = Q;
    rc(:,iter) = R;
    clearvars K Q R s2 s1;
    iter = iter +1;
    % sort all the pairs of this neuron, discarding any that have fewer spikes
%     [ccsort, ix] = sort(Xsim(isort(j),:) .* (nspk'>numel(s1)), 'descend');
%     ienu = find(ccsort<.5, 1) - 1; % find the first pair which has too low of a correlation

    % for all pairs above 0.5 correlation
%     for k = 1:ienu
%         s2 = rez.st3(rez.st3(:,2)==ix(k), 1)/ops.fs; % find the spikes of the pair
%         % compute cross-correlograms, refractoriness scores (Qi and rir), and normalization for these scores
%         [K, Qi, Q00, Q01, rir] = ccg(s1, s2, 40, dt);
%         Q = min(Qi/(max(Q00, Q01))); % normalize the central cross-correlogram bin by its shoulders OR by its mean firing rate
%         R = min(rir); % R is the estimated probability that any of the center bins are refractory, and kicks in when there are very few spikes
% 
%         if flag
%             if Q<.2 && R<.05 % if both refractory criteria are met
%                 i = ix(k);
%                 % now merge j into i and move on
%                 rez.st3(rez.st3(:,2)==isort(j),2) = i; % simply overwrite all the spikes of neuron j with i (i>j by construction)
%                 nspk(i) = nspk(i) + nspk(isort(j)); % update number of spikes for cluster i
%                 fprintf('merged %d into %d \n', isort(j), i)
%                 % YOU REALLY SHOULD MAKE SURE THE PC CHANNELS MATCH HERE
%                 break; % if a pair is found, we don't need to keep going (we'll revisit this cluster when we get to the merged cluster)
%             end
%         else
%           % sometimes we just want to get the refractory scores and CCG
%             rez.R_CCG(isort(j), ix(k)) = R;
%             rez.Q_CCG(isort(j), ix(k)) = Q;
% 
%             rez.K_CCG{isort(j), ix(k)} = K;
%             rez.K_CCG{ix(k), isort(j)} = K(end:-1:1); % the CCG is "antisymmetrical"
%         end
%     end
end
%toc;
% ccgout.xcorr = cc;
% ccgout.Q = qc;
% ccgout.R = rc;
ccgmetadata.lagamount = lagamount;
ccgmetadata.tbin = timebinsms;

% this part takes out autocorr from cc
if isequal(size(cc,1),size(cc,2))
    p = size(cc,2);     n = size(cc,3);
    ac = cc(bsxfun(@plus,[1:p+1:p*p]',[0:n-1]*p*p)); %#ok
    %from https://stackoverflow.com/questions/28402197/extract-diagonal-element-from-each-frontal-slice-of-tensor
else
    %batchlen = size(spiketimes,1);
    ac =zeros(size(cc,2),size(cc,3),'single');
    for kk = 1:size(cc,2), ac(kk,:) = squeeze(cc(batchlen(1)+kk-1,kk,:)); end
end

ccdata.crosscorr = cc;
ccdata.autocorr = ac;
ccdata.autocorrnopeak = ac;
ccdata.autocorrnopeak(ac==1) = 0;
ccdata.Q = qc;
ccdata.R = rc;
%ccdata.clusters = rawdata.clusters;
ccdata = struct2struct(ccdata,ccgmetadata);
ccdata.autocorrnopeak = ccdata.autocorrnopeak(:,ccdata.lagamount+1 :end );
ccdata.laghalfms = ccdata.lag(ccdata.lagamount+1:end)*1e6*ccdata.tbin;

end

