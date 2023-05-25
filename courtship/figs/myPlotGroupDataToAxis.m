function [ax, handles] = myPlotGroupDataToAxis(ax, ax_val, data, varargin)
% -------------------------------------------------------------------------
% function to hopefully streamline the plotting of grouped data by plotting
% just one data group to an axis. we'll always try to plot the raw data,
% and then give options for plotting either mean +/- se or median and IQR
%
% should probably figure out a better way to pass plot parrams, since the
% current way could lead to some headaches
%
% INPUTS:
%   - ax: handle of axis to plot to
%   - ax_val: value on axis to plot to (defaults to x axis value, but can
%       be y axis if "horzFlag" is true
%   - data: array of data to be plotted at ax_val
%   - varargnin: misc plot params
%       - plotColor: main color for data plot
%       - showMM: string to determine whether we should plot median and IQR
%           ('median'), mean +/- se ('mean'), or nothing ('none')
%       - jitterStrength: amount to jitter data by
%       - markerAlphaVal: transparency of raw data markers
%       - plotSpace: amount to shift mean plot by from data (if plotting)
%       - scatterSize: marker size for calls to scatter function
%       - markerSize: marker size for calls to plot function
%       - lineWidthThin: line width for things like mean errorbar
%       - lineWidthThick: line width for things like median line
%       - barWidth: width of bar for iqr patch
%       - barAlpha: transparency of bar for iqr patch
%       - barColor: color for iqr patch. Defaults to gray
%       - meanColor: color for mean +/- se errorbar plot (if using).
%           Defaults to black
%       - medianColor: color for line at median (if using). Defaults to
%           "plotColor" (defined above)
%       - horzFlag: boolean to flip to plotting data at a y axis value (as
%           opposed to x axis (default))
%       - capSize: cap size for errorbar plots
%       - fillFlag: for data scatter plot, fill points?
%       - scatterType: marker type for scatter plots
%       - meanType: marker type for mean errorbar plot
% -------------------------------------------------------------------------
%% parse inputs
if ~exist('ax','var') || isempty(ax)
    ax = gca ;
end
if ~exist('ax_val','var') || isempty(ax_val)
    ax_val = 1 ; % where to plot data
end
if ~exist('data','var') || isempty(data)
    data = nan ; % data to plot. if this is empty, just return
    fprintf('Error: no data provided. Quitting \n')
    return
end

% ----------------------------------
% defaults for plot params
defaults = {'showMM', 'median', ... % 'median' | 'mean'
    'plotColor', [0, 0, 0], ...  % default to black
    'jitterStrength', 0.25, ...
    'markerAlpha', 0.6, ...
    'plotSpace',0.33, ...
    'scatterSize', 8, ...
    'markerSize', 4, ...
    'lineWidthThin', 0.75, ...
    'lineWidthThick', 2.0 , ...
    'barWidth', 0.7, ...
    'barAlpha', 0.3 , ...
    'barColor', [0,0,0], ...
    'meanColor', [0, 0, 0], ...
    'medianColor', [], ...
    'horzFlag', false,...
    'capSize', 0, ...
    'fillFlag', false,...
    'scatterType', 'o',...
    'meanType', 'o'
    };

% ----------------------------------------------------------------------
% read variables or take default values
showMM          = get_var('showMM', 'defaults', defaults, varargin{:});
plotColor       = get_var('plotColor', 'defaults', defaults, varargin{:});
jitterStrength  = get_var('jitterStrength','defaults',defaults, varargin{:});
markerAlpha     = get_var('markerAlpha', 'defaults', defaults, varargin{:});
plotSpace       = get_var('plotSpace', 'defaults', defaults, varargin{:});
scatterSize     = get_var('scatterSize', 'defaults',defaults, varargin{:});
markerSize      = get_var('markerSize', 'defaults',defaults, varargin{:});
lineWidthThin   = get_var('lineWidthThin','defaults',defaults, varargin{:});
lineWidthThick  = get_var('lineWidthThick', 'defaults',defaults, varargin{:});
barWidth        = get_var('barWidth', 'defaults', defaults, varargin{:});
barAlpha        = get_var('barAlpha', 'defaults', defaults, varargin{:});
barColor        = get_var('barColor', 'defaults', defaults, varargin{:});
meanColor       = get_var('meanColor', 'defaults', defaults, varargin{:});
medianColor     = get_var('medianColor', 'defaults', defaults, varargin{:});
horzFlag        = get_var('horzFlag', 'defaults', defaults, varargin{:});
capSize         = get_var('capSize', 'defaults', defaults, varargin{:});
fillFlag        = get_var('fillFlag', 'defaults', defaults, varargin{:});
scatterType     = get_var('scatterType', 'defaults', defaults, varargin{:});
meanType        = get_var('meanType', 'defaults', defaults, varargin{:});

% -------------------------------------------
% axis properties (make sure hold is on)
holdFlag = ishold(ax) ;
hold(ax,'on') ;

% initialize array of handles
handles = gobjects(0) ;
% ---------------------------------------------------------
%% plot summary values for data

% plot either (median and iqr) or (mean +/- std) -- depends on showMM
switch showMM
    case 'median'
        % get median, 1st quartile, and 3rd quartile
        data_prctile = prctile(data, [25, 50, 75]) ;
        
        % get coordinates for box patch
        iqr_box_y = [data_prctile(1), data_prctile(3), data_prctile(3), ...
            data_prctile(1)] ;
        iqr_box_x = ax_val + [-barWidth, -barWidth, barWidth, ...
            barWidth]./2 ;
        
        % get coordinates for median line
        plot_x = ax_val + 0.5*barWidth.*[-1,1] ;
        plot_y = data_prctile(2).*[1,1] ;
        
        % switch x/y if we're plotting horizontally
        if horzFlag
            % iqr patch
            [iqr_box_x, iqr_box_y] = switchXY(iqr_box_x, iqr_box_y) ; 
           
            % median line
            [plot_x, plot_y] = switchXY(plot_x, plot_y) ; 
        end
        
        % draw iqr patch
        iqr_patch = patch(ax, iqr_box_x, iqr_box_y, barColor, ...
            'LineStyle', 'none') ;
        iqr_patch.FaceAlpha = barAlpha ;
        handles(end+1) = iqr_patch ; % add graphics object to handles array
        
        % draw line at median (need to check on color first though)
        if isempty(medianColor)
            medianColor = plotColor ; 
        end
        
        median_line = plot(ax, plot_x , plot_y, '-', ...
            'Color', medianColor,...
            'LineWidth',lineWidthThick) ;
        handles(end+1) = median_line ; % add graphics object to handles array
        
    case 'mean'
        % get data mean and standard DEVIATION
        data_mean = nanmean(data) ;
        data_std = nanstd(data) ;
        data_se = data_std / sqrt(length(data)) ;
        
        % define variables for x,y values to plot
        plot_x = ax_val+ plotSpace ;
        plot_y = data_mean ;
        yneg = data_se ; ypos = data_se ;
        xneg = [] ; xpos = [] ;
        
        % switch x/y for horizontal plot?
        if horzFlag
            % mean point
            [plot_x, plot_y] = switchXY(plot_x, plot_y) ; 
            
            % error bars
            [xneg, yneg] = switchXY(xneg, yneg) ; 
            [xpos, ypos] = switchXY(xpos, ypos) ; 
        end
        
        % plot mean +/- se to the right of where raw data will go
        mean_plot = errorbar(ax, plot_x, plot_y, yneg, ypos, xneg, xpos,...
            meanType,...
            'MarkerFaceColor',meanColor,...
            'Color',meanColor,...
            'MarkerSize',markerSize,...
            'CapSize',capSize,...
            'LineWidth',lineWidthThin) ;
        handles(end+1) = mean_plot ; % add graphics object to handles array
        
    case 'none'  
        % do nothing
    otherwise
        fprintf('Invalid selection for showMM \n')
        keyboard
end

% -------------------------------------------------------------------
%% plot raw data values
% get jitter for current data
jitter_curr = jitterStrength*(rand(size(data))-0.5) ;

% define x and y variables to plot
plot_x = ax_val*ones(size(data)) + jitter_curr ;
plot_y = data ;

% switch x,y if plotting horizontal?
if horzFlag
    [plot_x, plot_y] = switchXY(plot_x, plot_y) ; 
end

% fill in markers?
if fillFlag
   markerFaceColor = plotColor ; 
else
   markerFaceColor = 'none' ; 
end

% deal with the case where we have multiple scatter marker types (e.g. when
% we want to plot pitch up vs down with '^' vs 'v')
if iscell(scatterType)
    % get unique scatter marker types
    [uniqueTypes, ~, typeIdx] = unique(scatterType) ; 
    
    % loop over marker types and plot for each
    for k = 1:length(uniqueTypes)
        % current marker type
        scatterTypeCurr = uniqueTypes{k} ; 
        
        % index of data values with current marker type
       idx = (typeIdx == k) ; 
       
       % make plot for current markers
       data_scatter = scatter(ax, plot_x(idx), plot_y(idx), scatterSize, ...
           plotColor, scatterTypeCurr,...
           'MarkerEdgeColor', plotColor, ...
           'MarkerEdgeAlpha', markerAlpha,...
           'MarkerFaceColor', markerFaceColor, ...
           'MarkerFaceAlpha', markerAlpha) ;
       handles(end+1) = data_scatter ; % add graphics object to handles array
    end
else
    % otherwise just call scatter once to plot raw data 
    data_scatter = scatter(ax, plot_x, plot_y, scatterSize, ...
        plotColor, scatterType,...
        'MarkerEdgeColor', plotColor, ...
        'MarkerEdgeAlpha', markerAlpha,...
        'MarkerFaceColor', markerFaceColor, ...
        'MarkerFaceAlpha', markerAlpha) ;
    handles(end+1) = data_scatter ; % add graphics object to handles array
end

% ----------------------------------------------------
%% return axis properties to previous values
if ~holdFlag
    hold(ax,'off') ;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HELPER FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------
%% function to switch values of two variables
function [x_out, y_out] = switchXY(x_in, y_in)
x_out = y_in ;
y_out = x_in ;
end