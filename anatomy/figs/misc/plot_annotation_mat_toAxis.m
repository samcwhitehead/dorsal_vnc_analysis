% -------------------------------------------------------------------------
% function to plot annotation matrix to an extant axis
% -------------------------------------------------------------------------
function ax = plot_annotation_mat_toAxis(fig, ax, annotationVals, ...
    sort_ind, y_tick_labels, x_tick_labels, flipYTickLabelFlag, color_mat,...
    ylim, xlim, fontSize, fontSizeSmall, fontName, gridFlag)
% ---------------------------------------------
%% inputs and params
% parent figure window
if ~exist('fig','var') || isempty(fig)
   fig = gcf ;  
end
% axis to plot to
if ~exist('ax','var') || isempty(ax)
   ax = gca ;  
end
% % reshape annotationVals if need be
% if size(annotationVals,1) < size(annotationVals,2)
%     annotationVals = annotationVals' ;
% end
% sorting for annotation values (sort neurons)
if ~exist('sort_ind','var') || isempty(sort_ind) 
   sort_ind = 1:size(annotationVals,1) ; 
end
% target labels
if ~exist('y_tick_labels','var') || isempty(y_tick_labels)
   y_lick_labels = 1:size(annotationVals,2) ;  
end
% neuron labels
if ~exist('x_tick_labels','var') || isempty(x_tick_labels)
   x_lick_labels = 1:size(annotationVals,1) ;  
end
% whether or not to flip y axis labels (so that whole figure can be
% rotated 90 degrees)
if ~exist('flipYTickLabelFlag','var') || isempty(flipYTickLabelFlag)
   flipYTickLabelFlag = true ;  
end
% colors to define input, output, and both
if ~exist('color_mat','var') || isempty(color_mat)
    color_mat = [1.0 , 1.0, 1.0 ;  [33,102,172]/255 ; [178,24,43]/255 ; ...
        [118,42,131]/255] ;
end
% axis y limits
if ~exist('ylim','var') || isempty(ylim)
    ylim = [0.5, length(y_tick_labels) + 0.5] ;
end
% axis x limits
if ~exist('xlim','var') || isempty(xlim)
    xlim = [0.5, length(x_tick_labels) + 0.5] ; 
end
% general font size
if ~exist('fontSize','var') || isempty(fontSize)
    fontSize = 6 ; 
end
% small font size to use for dense labels
if ~exist('fontSizeSmall','var') || isempty(fontSizeSmall)
    fontSizeSmall = 4.5 ;
end
% which font to use
if ~exist('fontName','var') || isempty(fontName)
    fontName = 'arial' ; 
end
% include grid for matrix plot?
if ~exist('gridFlag','var') || isempty(gridFlag)
   gridFlag = true;  
end

% check properties of current axis
set(fig,'CurrentAxes',ax)
hold(ax,'on')

% line width for grid
lw_grid = 0.5 ;

% ---------------------------------------------------
%% make sure neuron labels are formatted correctly
for n = 1:length(x_tick_labels)
   x_tick_labels{n} = strrep(x_tick_labels{n},'_','-') ; 
end
% -----------------------------------------------
%% generate plot
% sort annotation values if we're plotting cluster results
annotationVals_sort = annotationVals(sort_ind,:) ;

% adjust colormat to only include values that we have present in
% annotationVals
annotationVals_unique = unique(annotationVals(:)) ; 
if ~all(round(annotationVals_unique) == annotationVals_unique)
    % account for 1-start indexing and fractional values
    color_ind = 2*annotationVals_unique + 1 ; 
else
    color_ind = annotationVals_unique + 1 ; % account for 1-start indexing
end
%color_mat_new = ones(size(color_mat)) ; 
%color_mat_new(color_ind,:)  = color_mat(color_ind,:) ; 
color_mat_new = color_mat(1:max(color_ind),:) ; 

% --------------------
% plot image to axis
imagesc(ax, flipud(annotationVals_sort'))
colormap(ax, color_mat_new)

% -------------------------
% plot grid lines
if gridFlag
    for i = 0:(size(annotationVals_sort,1))
        plot(ax, [i + 0.5, i + 0.5], [0.5, size(annotationVals_sort,2)+0.5], ...
            'k-', 'LineWidth',lw_grid)
    end
    for j = 0:size(annotationVals_sort,2)
        plot(ax, [0.5,size(annotationVals_sort,1)+0.5], [j + 0.5, j + 0.5 ], ...
            'k-', 'LineWidth',lw_grid)
    end
end
% ---------------------------
% set axis properties
%axis equal
set(ax,'xlim',xlim,'ylim',ylim)

% axis labels
set(ax,'fontName',fontName)
ax.XAxis.FontSize = fontSizeSmall ;
ax.YAxis.FontSize = fontSize ;

set(ax,'XTick',1:size(annotationVals_sort,1), ...
    'XTickLabel',x_tick_labels(sort_ind),...
    'XTickLabelRotation',90)
set(ax,'YTick',1:size(annotationVals_sort,2), ...
    'YTickLabel',fliplr(y_tick_labels))
if flipYTickLabelFlag
    set(ax,'YTickLabelRotation',180) ;
end
pause(0.1)
try
    ax.XRuler.Axle.LineStyle = 'none' ; 
catch
    fprintf('No x axis ruler... \n')
end
try
    ax.YRuler.Axle.LineStyle = 'none' ; 
catch
    fprintf('No y axis ruler... \n')
end
    
set(ax,'TickLength',[0 0])
end