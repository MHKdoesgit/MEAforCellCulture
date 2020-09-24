

function varargout = rasterPlotter(spikes,varargin)
%
%%% rasterPlotter %%%
%
%
% This function make a rasters plot with lines like classical papers.
%
% ===============================Inputs====================================
%
%   spikes : spikes responses to stimulus.
%   plot_color : plot color (optional).
%   plot_linwidth : thickness of each spike line.
%   spikes_height : height of each spike line.
%
%================================Output====================================
%
%   plot : rasters plot of the spikes.
%
%   caution: This function is drawing individual line by looping over all
%   the spikes. Only use it for aesthetic since it is very slow.
%
% written by Mohammad, 29.10.2015
% update to superfast version with lineplot on 31.03.2016 (best
% procastination ever! ;)).
% improvment on the speed on line plot by using single line instead of
% multiple lines inspired by
% http://undocumentedmatlab.com/blog/some-performance-tuning-tips on
% 11.04.2016


if (nargin > 1 && isnumeric(varargin{1})), spikes_height = varargin{1};  else, spikes_height = 0.4;  end
if nargin > 2, plot_color = varargin{2};  else, plot_color = [0 0.7461 1]; end
if nargin > 3, plot_linwidth = varargin{3};  else, plot_linwidth = 1;  end

% if size(spikes,1) > size(spikes,2)
%      spikes = transpose(spikes);
% end;

% if the input is binary with zeros and ones find the time of spikes
if isempty(spikes(spikes~=0 & spikes~=1))
    spikes = getSpiketimesfromLogical(spikes);
    if nargout == 2
        varargout{2} = spikes;
    end
end

% error checking
if not(isnumeric (spikes_height)),  spikes_height = 0.4; end
if length(spikes_height) ~=1, spikes_height = 0.4; end
% making top and buttom point of the spikes
spktp = repmat(1:size(spikes,2),size(spikes,1),1) + spikes_height;
spkbp = repmat(1:size(spikes,2),size(spikes,1),1) - spikes_height;

% in order to reduce the number of displayed graphic handles, we can unify the separate line
% segments into a single line that has NaN (or Inf) values interspaced between the segments.
spkoneline = [spikes(:),spikes(:),nan(numel(spikes(:)),1)]';
tboneline = [spkbp(:),spktp(:), nan(numel(spktp(:)),1)]';

try   % small trick to make it more flexiable to user inputs
    rasplt = line(spkoneline(:),tboneline(:),'color',plot_color,'linewidth',plot_linwidth);
catch ME 
    disp(ME.message)
    rasplt = plot(spkoneline(:),tboneline(:),varargin{:});
end

hold off;
if nargout == 1
varargout{1} = rasplt;
end

% old and slow method
% rasplt = line([spikes(:),spikes(:)]', [spkbp(:),spktp(:)]','color',plot_color,'linewidth',plot_linwidth);
% rasplt = plot([spikes(:),spikes(:)]', [spkbp(:),spktp(:)]',varargin{:});

%%% old way of plotting, very slow

% trials = size(spikes,1);
% hold on;
% for j= 1:trials
%     spiketime = spikes(j,~isnan(spikes(j,:)));
%     spkposX = [spiketime;spiketime];
%     spkposY = [((1:length(spiketime))-spikes_height); ((1:length(spiketime))+spikes_height)];
%     line(spkposX,spkposY,'color',plot_color,'linewidth',plot_linwidth);
%
% %     pl = line([spiketime(1),spiketime(1)],[j-spikes_height, j+spikes_height],'color',plot_color,'linewidth',plot_linwidth);
% %     %hold on;
% %     drawnow;
% %     for i = 2:length(spiketime),
% %         set(pl,'Xdata',[spiketime(i),spiketime(i)], 'Ydata',[j-spikes_height, j+spikes_height]);
% %         drawnow;
% %     end;
%     clear spiketime spkposX spkposY;
% end;
% hold off;

end

%--------------------------------------------------------------------------------------------------%
%----------                               sub-functions                                  ----------%
%--------------------------------------------------------------------------------------------------%

function spktimes = getSpiketimesfromLogical(spk)
spkidx = cell(1,size(spk,1));
for i = 1:size(spk,1)
    spkidx{i} = find(spk(i,:) == 1);
end
spktimes = CelltoMatUE(spkidx);

end