

function b = barwitherror(d,colset,lw, varargin)
%
%%% thesis.barwitherror %%%
%
%
% This function plots the input data in form of nice bars with the error
% bar as SEM on top of each bar.
%
%================================Inputs====================================
%
%   data : input data for multiple boxplot in cell format.
%   colset : color set for each boxplot.
%   lw : line width for the error bars.
%
%================================Output====================================
%
%   No output : since this a just a plotting function there is no output
%   but in case of request the handels of the plot are returened.
%
% written by Mohammad, 09.08.2017, moved to separated function on 24.08.2017
% added the flag for drawing the data points on 11.09.2017.
% updated to new MATLAB version on 11.01.2019.


if nargin > 3, drawdataFlag = varargin{1}; else, drawdataFlag = false; end
if nargin > 4, capsize = varargin{2}; else, capsize = 8; end
if nargin > 5, alph = varargin{3}; else, alph = 0.8; end

ymean = nanmean(CelltoMatUE(d));
ysem = nansemSuite(CelltoMatUE(d));
x = 1:length(ymean);
b = zeros(length(x),1);

if drawdataFlag     % flag for drawing the data points next to each bar
    pos = 1:0.5:1.5*size(d,2)+0.5;
    boxpos = pos(2:3:end);
    sprdpos = pos(1:3:end)-0.15;
    
    sprdcol = mat2cell(colset,ones(size(d,2),1));
    s = plotSpread(CelltoMatUE(d),'spreadwidth',0.45,'distributionMarkers','o','distributionColors',...
        sprdcol,'xValues',sprdpos);
    x = boxpos;
end

for ii = 1:length(x)
    
    b = bar(x(ii),ymean(ii));
    if not(ishold),     hold on;    end
    b.FaceColor = colset(ii,:); b.EdgeColor = 'none';
    b.FaceAlpha = alph;
    %set(b(ii),'facecolor',colset(ii,:),'edgecolor','none');
    %set(get(b,'Children'),'FaceAlpha',0.3)
    
    %c = get(b(ii),'Children');
   %  xdata = mean(b.XData);
   %  tempYData  = b.YData;
   %  ydata = mean(tempYData);%mean(tempYData(2:3,:))';
    
    %ha = errorbar(xdata,ydata,ystd(ii),'color',colset(ii,:),'linewidth',lw);
    ha = errorbar(b.XData,mean(b.YData),ysem(ii),'color',colset(ii,:),'linewidth',lw);
    
    % this is from matlab answers to change the capsize of the error bars.
    %hb = get(ha,'children');
    %Xdata = ha.XData;%get(hb(2),'Xdata');
    %temp = 4:3:length(Xdata);
    %temp(3:3:end) = [];
    % xleft and xright contain the indices of the left and right
    %  endpoints of the horizontal lines
    %xleft = temp; xright = temp+1;
    % Increase line length by 0.2 units
    %Xdata(xleft) = Xdata(xleft) - .1;
    %Xdata(xright) = Xdata(xright) + .1;
    %set(hb(2),'Xdata',Xdata);
    
    ha.CapSize = capsize;
    
    if drawdataFlag
    set(s{1}(ii),'markerfacecolor', colset(ii,:), 'markeredgecolor',abs(colset(ii,:)-0.2),'markersize',3);
    end
    
    % clearvars c e xdata ydata tempYData ha hb Xdata temp xleft xright;
end

if drawdataFlag
    xlim([0 max(pos)+0.5]);
    set(gca,'xtick',(sprdpos + boxpos)/2,'xticklabel',1:length(x),'tickdir','out');
else
    set(gca,'xtick',x,'tickdir','out');
end
% if nargout > 1
%     varargout{1} = b;
% end;

end

