% -------------------------------------------------------------------------
% function to plot overlap results in line form
% -------------------------------------------------------------------------
function h_overlapLine = plot_overlap_data_lines(overlap_data, labels_x, ...
    labels_y, plotColor, axisToCompress, overlap_lim)
% -------------------------
%% params and inputs
if ~exist('labels_x','var') || isempty(labels_x) 
    labels_x = [] ; 
end
if ~exist('labels_y','var') || isempty(labels_y) 
    labels_y = [] ; 
end
if ~exist('plotColor','var') || isempty(plotColor)
    plotColor = [0, 0, 0] ;
end
if ~exist('axisToCompress','var') || isempty(axisToCompress)
    axisToCompress = 'X' ;
end
if ~exist('overlap_lim','var') || isempty(overlap_lim)
    overlap_lim = [0, 0.15] ;
end
% -------------------
% plot preferences

lineAlpha = 0.35 ; 
lineWidthThin = 0.5 ; 
lineWidthThick = 1.5 ; 

axisFontSize = 6 ; % needs to be tiny for labels
axisFontSizeSmall = 4.5 ; % needs to be tiny for labels

plotColorMean = [0, 0, 0] ; 

% number of stacks compared
[N_stacks1, N_stacks2] = size(overlap_data) ;

% ---------------------------------------------------------------
%% get overlap values (different for table vs mat)
if istable(overlap_data)
    overlap_vals = overlap_data{:,:} ;
else
    overlap_vals = overlap_data ;
end

% -----------------------------------------------------------
%% make figure
h_overlapLine = figure('PaperPositionMode','auto') ;
hold on
switch axisToCompress
    case {'x','X'}
        % lines for individual neurons
        for ind = 1:N_stacks2
           plot(overlap_vals(:,ind), 1:N_stacks1, '-', ...
               'Color', [plotColor, lineAlpha], 'LineWidth', lineWidthThin) 
        end
        % line for median
        mean_overlap = nanmean(overlap_vals,2) ; 
        plot(mean_overlap, 1:N_stacks1,  '-','Color', plotColorMean, ...
            'LineWidth', lineWidthThick) 
        
        % axis properties
        ax = gca ;
        
        set(ax, 'xlim', overlap_lim)
        set(ax, 'ylim', [1, N_stacks1])
        set(ax,'YTick',1:N_stacks1, 'YTickLabel',labels_y)
        set(ax,'ydir','reverse')
        
        % shorten/remove ticks
        ax.XAxis.TickLength = [0.01, 0.025] ; 
        ax.YAxis.TickLength = [0.0, 0.0] ; 
        
    case {'y','Y'}
        % lines for individual neurons
        for ind = 1:N_stacks1
           plot(1:N_stacks2, overlap_vals(ind,:), '-', ...
               'Color', [plotColor, lineAlpha], 'LineWidth', lineWidthThin) 
        end
        % line for median
        mean_overlap = nanmean(overlap_vals,1) ; 
        plot(1:N_stacks2, mean_overlap, '-','Color', plotColorMean, ...
            'LineWidth', lineWidthThick) 
        
        % axis properties
        ax = gca ;
        set(ax, 'xlim', [1, N_stacks2])
        set(ax,'XTick',1:N_stacks2, 'XTickLabel',labels_x)
        set(ax, 'ylim', overlap_lim)
        
         % shorten/remove ticks
        ax.XAxis.TickLength = [0.0, 0.0] ; 
        ax.YAxis.TickLength = [0.01, 0.025] ; 
    otherwise
        fprintf('Invalid axis choice: %s \n', axisToCompress)
        keyboard
end

% general axis properties
set(ax,'fontName','arial','fontSize',axisFontSize,...
    'TickLabelInterpreter','none')

if (N_stacks2 > 100)
    ax.XAxis.FontSize = axisFontSizeSmall ;
end
if (N_stacks1 > 100)
    ax.YAxis.FontSize = axisFontSizeSmall ;
end



end 

