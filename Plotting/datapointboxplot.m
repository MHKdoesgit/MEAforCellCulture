

function s = datapointboxplot(data, varargin)
%
%%% datapointboxplot %%%
%
%
% This function plots the input data in form of nice boxplots with the real
% distribution of the data next to each boxplot.
%
%================================Inputs====================================
%
%   data : input data for multiple boxplot in cell format.
%   colset : color set for each boxplot.
%   dpointsflag : flag to draw points on top of each boxplot.
%
%================================Output====================================
%
%   No output : since this a just a plotting function there is no output
%   but in case of request the handels of the plot are returened.
%
% written by Mohammad, 09.08.2017, moved to separated function on
% 24.08.2017.
% major update by adding more features and optional input plus transparent
% points and boxes, on 08.02.2019

if iscell(data),    data = CelltoMatUE(data);       end
if size(data,1)<size(data,2), data = transpose(data); end
%defaults:____________________
def.colset = hsv(size(data,2)); %set of colors for each box.
def.dpointsflag = false; % flag for showing data on top of each box.
def.boxwidth = 0.5;    % box width.
def.spreadwidth = def.boxwidth/2;   % width of data spread next to each box.
def.boxspreadgap = 0; % the relative gap between each box and spread point next to it.
def.boxlw = 0.5; % line width of each boxes
def.boxmedianlw= def.boxlw;     % width of the median line in the middle of each box.
def.scatterlw = 0.001; % Marker line width for scatter data.
def.boxalpha = 0.2;  % box transparency.
def.scattersize = 15;   % scatter Marker size.
def.scfacealpha = 0.4;   % scatter face alpha.
def.scedgealpha = 0.5;  % scatter edge alpha.

% parsing inputs: _____________________
pnames = {'Colors','ShowDatapoints','BoxWidth','SpreadWidth','bxspgap','BoxLinewidth','BoxMedianLinewidth','ScatterLinewidth',...
    'BoxAlpha','ScatterMarkerSize','ScatterFaceAlpha','ScatterEdgeAlpha'};

[para.colset, para.dpointsflag, para.boxwidth, para.spreadwidth, para.bxspgap ,para.boxlw, para.boxmedianlw, para.scatterlw,para.boxalpha, ...
    para.scattersize, para.scfacealpha, para.scedgealpha] = internal.stats.parseArgs(lower(pnames), struct2cell(def)', varargin{:});

pos = (1:3*size(data,2)+1)/2;%(1:0.5:1.5*size(data,2)+0.5);
boxpos = pos(1:3:end);  boxpos = boxpos(1:size(data,2));
sprdpos = pos(2:3:end)-para.bxspgap;   % sprdpos = sprdpos(1:length(boxpos));

% boxpos = pos(2:3:end);  
% sprdpos = pos(1:3:end)-para.bxspgap;    sprdpos = sprdpos(1:length(boxpos));

stemp = plotSpread(data,'spreadwidth',para.spreadwidth,'xValues',sprdpos);
s.xData = get(stemp{:,1},'XData');    s.yData = get(stemp{:,1},'YData');
delete(findobj(stemp{3}, 'type', 'line'));
boxplot(data,'positions',boxpos,'boxstyle', 'outline','widths', para.boxwidth,'symbol',...
    'w.','color',para.colset, 'OutlierSize',0.1);
hold on;
% sprdcol = Piotr_Toolbox.matlab.mat2cell2(colset(1:size(data,2),:),[size(data,2),1]);
% s = plotSpread(data,'spreadwidth',boxwidth/2,'distributionMarkers','o','distributionColors',...
%     sprdcol,'xValues',sprdpos);
set(gca,'xtick',(sprdpos + boxpos)/2-(para.bxspgap/2),'tickdir','out');
xlim([sprdpos(1)-sprdpos(1), boxpos(end)+sprdpos(1)])

boxes = findall(gcf, 'tag', 'Box');
bxoutliers = flipud(findall(gcf, 'tag', 'Outliers'));
bxwiskersUP = findall(gcf, 'tag', 'Upper Whisker');
bxwiskersDN = findall(gcf, 'tag', 'Lower Whisker');
bxAdjUP = findall(gcf, 'tag', 'Upper Adjacent Value');
bxAdjDN = findall(gcf, 'tag', 'Lower Adjacent Value');
bxmedians = findall(gcf, 'tag', 'Median');

for ii = 1:size(data,2)
    patch(get(boxes(1+size(data,2)-ii),'xData'),get(boxes(1+size(data,2)-ii),'yData'),1,...
        'Facecolor',para.colset(ii,:),'Edgecolor','none','Facealpha',para.boxalpha);
    set(boxes(ii),'LineWidth', para.boxlw);
    set(bxwiskersUP(ii),'LineWidth', para.boxlw,'linestyle','-');
    set(bxwiskersDN(ii),'LineWidth', para.boxlw,'linestyle','-');
    set(bxAdjUP(ii),'LineWidth', para.boxlw,'linestyle','-');
    set(bxAdjDN(ii),'LineWidth', para.boxlw,'linestyle','-');
    set(bxmedians(ii),'LineWidth', para.boxmedianlw);
    %set(bxoutliers(ii),'markerfacecolor', colset(1+size(d,2)-ii,:), 'markeredgecolor','none');
    set(bxoutliers(ii),'markerfacecolor', abs(para.colset(ii,:)-0.2), 'markeredgecolor','none');
    %  set(stemp{1}(ii),'markerfacecolor', colset(ii,:), 'markeredgecolor',abs(colset(ii,:)-0.2),'markersize',3);
    if not(iscell(s.xData))
        sx = s.xData;   sy = s.yData;
    else
        sx = s.xData{ii};   sy = s.yData{ii};
    end
    
    scatter(sx,sy,para.scattersize,'filled','MarkerFaceColor',para.colset(ii,:),...
        'MarkerFaceAlpha',para.scfacealpha,'markeredgecolor',abs(para.colset(ii,:)-0.2),'MarkerEdgeAlpha',...
        para.scedgealpha,'LineWidth',para.scatterlw);
    
    if para.dpointsflag
        hold on;
        plot(boxpos(ii),data(:,ii),'o','markersize',sqrt(para.scattersize),'markerfacecolor', ...
            para.colset(ii,:),'markeredgecolor', 'none')
    end
end

end




