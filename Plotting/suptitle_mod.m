

function hout=suptitle_mod(currentfig, str,enlarge)
%
%SUPTITLE Puts a title above all subplots.
%	SUPTITLE('text') adds text to the top of the figure
%	above all subplots (a "super title"). Use this function
%	after all subplot commands.

% Drea Thomas 6/15/95 drea@mathworks.com
% Adjusted by Mohammad to get the figure and font size

% Warning: If the figure or axis units are non-default, this
% will break.


if not(ishandle(currentfig)),
    str = currentfig;
    currentfig = gcf;
end;

if nargin < 3, enlarge = 2; end;

% Parameters used to position the supertitle.

% Amount of the figure window devoted to subplots
plotregion = .92;

% Y position of title in normalized coordinates
titleypos  = .95;

% Fontsize for supertitle
fs = get(currentfig,'defaultaxesfontsize')+ enlarge;

% Fudge factor to adjust y spacing between subplots
fudge=1;

%haold = currentfig;
figunits = get(currentfig,'units');

% Get the (approximate) difference between full height (plot + title
% + xlabel) and bounding rectangle.

if (~strcmp(figunits,'pixels')),
    set(currentfig,'units','pixels');
    pos = get(currentfig,'position');
    set(currentfig,'units',figunits);
else
    pos = get(currentfig,'position');
end
ff = (fs-4)*1.27*5/pos(4)*fudge;

% The 5 here reflects about 3 characters of height below
% an axis and 2 above. 1.27 is pixels per point.

% Determine the bounding rectange for all the plots

% h = findobj('Type','axes');

% findobj is a 4.2 thing.. if you don't have 4.2 comment out
% the next line and uncomment the following block.

h = findobj(currentfig,'Type','axes');  % Change suggested by Stacy J. Hills

max_y=0;
min_y=1;

oldtitle =0;
for i=1:length(h),
    if (~strcmp(get(h(i),'Tag'),'suptitle')),
        pos=get(h(i),'pos');
        if (pos(2) < min_y), min_y=pos(2)-ff/5*3;end;
        if (pos(4)+pos(2) > max_y), max_y=pos(4)+pos(2)+ff/5*2;end;
    else
        oldtitle = h(i);
    end
end

if max_y > plotregion,
    scale = (plotregion-min_y)/(max_y-min_y);
    for i=1:length(h),
        pos = get(h(i),'position');
        pos(2) = (pos(2)-min_y)*scale+min_y;
        pos(4) = pos(4)*scale-(1-scale)*ff/5*3;
        set(h(i),'position',pos);
    end
end

np = get(currentfig,'nextplot');
set(currentfig,'nextplot','add');
if (oldtitle),
    delete(oldtitle);
end
ha = axes('pos',[0 1 1 1],'visible','off','Tag','suptitle');
ht = text(.5,titleypos-1,str);
set(ht,'horizontalalignment','center','fontsize',fs);
set(currentfig,'nextplot',np);
%axes(haold);
if nargout,
    hout=ht;
end