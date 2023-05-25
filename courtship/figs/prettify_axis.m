% -------------------------------------------------------------------------
% function to (ideally) automate some figure formatting stuff so that i
% don't have to mess around in inkscape as much
% -------------------------------------------------------------------------
function ax = prettify_axis(ax, offsetFlag, symXFlag, symYFlag, ...
    sparseXTickFlag, sparseYTickFlag, labelDistShrink )
% ---------------------------------
%% inputs and params
if ~exist('ax','var') || isempty(ax)
    ax = gca ; 
end
if ~exist('offsetFlag','var') || isempty(offsetFlag)
    offsetFlag = true ; 
end
if ~exist('symYFlag','var') || isempty(symYFlag)
    symYFlag = false ; 
end
if ~exist('symXFlag','var') || isempty(symXFlag)
    symXFlag = true ; 
end
if ~exist('sparseXTickFlag','var') || isempty(sparseXTickFlag)
    sparseXTickFlag = true ; 
end
if ~exist('sparseYTickFlag','var') || isempty(sparseYTickFlag)
    sparseYTickFlag = true ; 
end
if ~exist('labelDistShrink','var') || isempty(labelDistShrink)
    labelDistShrink = 0.85 ; % amount to multiply distance between axis label and line
end
% import general plot preferences
b1_paper_plot_pref ; 

% other general params:

% ------------------------------------------
%% first set axis font and line properties
set(ax,'LineWidth', axisLineWidth,...
    'fontName', fontName, ...
    'fontSize',axisFontSize,...
    'fontWeight','normal')

% if there are axis labels, adjust them too
xlab = get(ax,'xlabel') ;
set(xlab,'fontName', fontName, ...
        'fontSize',labelFontSize,...
        'fontWeight','normal')


ylab = get(ax,'ylabel') ;
set(ylab,'fontName', fontName, ...
    'fontSize',labelFontSize,...
    'fontWeight','normal')

% plus title
tit = get(ax,'title') ;
set(tit,'fontName', fontName, ...
    'fontSize',titleFontSize,...
    'fontWeight','normal')

% to self: title as well?
% --------------------------------------------------
%% now deal with limits and tick marks
% first tighten but add a bit of white space
%ax = axisNotSoTight(ax) ; 

% if we want the axes symmetric, make it so:
if symXFlag
   xmax = max(abs(ax.XLim)) ; 
   set(ax,'xlim',xmax*[-1,1]) ; 
end

if symYFlag
   ymax = max(abs(ax.YLim)) ; 
   set(ax,'ylim',ymax*[-1,1]) ; 
end

% set axis tick direction
ax.TickDir = 'out' ; 

% this is where it might get a little goofy. now want to perform axis
% offset (i.e. separate x and y axes at ends) if selected. then we'll trim
% the tick marks to make sure they don't look weird or are too dense
if offsetFlag
    offsetAxes(ax) ;
end

% sparsen x ticks
if sparseXTickFlag 
    x_ticks = get(ax,'XTick') ; 
    % x ticks 
    if symXFlag 
       x_ticks_new = [x_ticks(1), 0, x_ticks(end)] ; 
    elseif (0 > x_ticks(1)) && (0 < x_ticks(end))
        if (x_ticks(end) >= 3*abs(x_ticks(1)))
            x_ticks_new = [0, x_ticks(end)] ; 
        else
            tick_min_x = min(abs([x_ticks(1), x_ticks(end)])) ;  
            x_ticks_new = [-1*tick_min_x, 0, tick_min_x] ; 
        end
    else
        x_ticks_new = [x_ticks(1), x_ticks(end)] ; 
    end
    
    
    % set new tick labels
    set(ax,'XTick',x_ticks_new) ; 
    
end

% sparsen y ticks
if sparseYTickFlag
    y_ticks = get(ax,'YTick') ;
    
    % y ticks
    if symYFlag
        tick_min_y = min(abs([y_ticks(1), y_ticks(end)])) ;
        y_ticks_new = [-1*tick_min_y, 0, tick_min_y] ;
        
    elseif (0 > y_ticks(1)) && (0 < y_ticks(end))
        if (y_ticks(end) >= 3*abs(y_ticks(1)))
            y_ticks_new = [0, y_ticks(end)] ;
        else
%             y_ticks_new = [y_ticks(1), 0, y_ticks(end)] ;
            y_ticks_new = [y_ticks(1), y_ticks(end)] ;
        end
    else
        y_ticks_new = [y_ticks(1), y_ticks(end)] ;
    end
    
    set(ax,'YTick',y_ticks_new) ;
end

% -----------------------------------
%% alter label to axis distance
% check if we need different values for x and y labels
if (length(labelDistShrink) > 1)
    labelDistShrinkX = labelDistShrink(1) ;
    labelDistShrinkY = labelDistShrink(2) ;
else
    labelDistShrinkX = labelDistShrink ;
    labelDistShrinkY = labelDistShrink ;
end
if ~strcmp(xlab.String,'')
    %move x labels a teensy bit closer to axis line
    set(xlab,'Units','normalized')
    xlab_pos = get(xlab,'Position') ;
    xlab_pos(2) = labelDistShrinkX*xlab_pos(2) ;
    set(xlab,'Position',xlab_pos)
end

if ~strcmp(ylab.String,'')
    % move y labels a teensy bit closer to axis line
    set(ylab,'Units','normalized')
    ylab_pos = get(ylab,'Position') ;
    ylab_pos(1) = labelDistShrinkY*ylab_pos(1) ;
    set(ylab,'Position',ylab_pos)
end

end