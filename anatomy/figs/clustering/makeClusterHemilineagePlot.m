% -------------------------------------------------------------------------
% script to generate grid of pie charts showing hemilineage breakdown of
% each IN cluster. going to borrow code from the digraph plots so that the
% two panels look similar
% -------------------------------------------------------------------------
% path info 
[mfilePath, ~, ~] = fileparts(mfilename('fullpath')) ; 
figDirectory = fileparts(mfilePath) ;
parentDirectory = fileparts(figDirectory) ; 
rootPath = fullfile(parentDirectory, 'data', 'clustering') ;
savePath = fullfile(mfilePath,'output') ;

% file with clustering results
dataFilename =  'conn_clust_struct_IN.mat' ; 

% save results?
saveFlag = false ; % for plots

% plot results as bars or pie charts?
plotType = 'bar' ; % 'pie' | 'bar'

% how to color plots ('hemilineage' | 'behavior' | 'neurotransmitter' |
% 'cluster')
colorMode = 'neurotransmitter' ;  % 'neurotransmitter' ;  %'behavior' ; % 'hemilineage' ;

% -------------------------
% general plot params
figUnits = 'inches' ;
legendFlag = ~strcmp(plotType,'bar') ; % haven't added this to bar plot option yet

% pie chart params
useSubTypeFlag = false ;
titleFlag = true ;
scaleChartFlag = true ;
sortWedgeFlag = true ;

% bar plot params
sortStackFlag = true ;
fullLabelFlag = true ;
horizontalFlag = false ;

% -----------------------------------------------------------
%% change save path depending on plot type
switch plotType
    case 'pie'
        savePath = fullfile(savePath,'pie_charts') ; 
    case 'bar'
        savePath = fullfile(savePath,'bar_plots') ; 
    otherwise
        fprintf('Error: invalid plot type \n')
        keyboard
end

% make sure save directory exists
if ~exist(savePath,'dir')
   mkdir(savePath) 
end
% ------------------------------------------------------------------------
%% load previous clustering results
connPath = fullfile(rootPath, dataFilename) ;
clustDir = dir(connPath) ; 
if length(clustDir) == 1
    conn_clust_struct = importdata(fullfile(clustDir(1).folder, ...
        clustDir(1).name)) ;
else
    fprintf('No clustering results found -- quitting \n')
    return
end

clustNums = unique([conn_clust_struct(:).cluster_idx]) ;
N_clusts = length(clustNums) ;

%------------------------------------------------------------------
%% switch plot params based on plot type
switch plotType
    case 'pie'
        % subplot tiling dimensions
        n_cols = 4 ;
        n_rows = ceil(N_clusts/n_cols) ;
        
        
        % figure outer position
        figPosition = [0.5, 2.3229, 6.5, 1.6875*n_rows + 0.45 ] ;   %[23.7396    2.3229    4.9    5.6667] ;   % [24.0000    2.7604   5.0313    5.0938]
        
        % axis position (pies):
        gap = [0.04, 0.01] ; %[0.03, 0.05] ; %[0.1 0.10] ;
        marg_h = [0.03, 0.03] ; %[0.08 0.08] ;
        marg_w = [0.03, 0.03] ;
        
        % x,y axis limits for each pie chart
        xlim = 1.4*[-1, 1] ;
        ylim = 1.2*[-1, 1] ;
        
    case 'bar'
        % figure outer position
        if horizontalFlag
            figPosition = [0.5, 2.3229, 3.65, 6.5] ;
            % axis position (bar plot). NB: this doesn't matter a ton, since
            % all clusters are on the same axis
            gap = [0.01, 0.01] ; %[0.03, 0.05] ; %[0.1 0.10] ;
            marg_h = [0.05, 0.03] ;
            marg_w = [0.12, 0.01] ; 
        else
            figPosition = [0.5, 2.3229, 6.5, 3.0] ;
            % axis position (bar plot). NB: this doesn't matter a ton, since
            % all clusters are on the same axis
            gap = [0.01, 0.01] ; %[0.03, 0.05] ; %[0.1 0.10] ;
            marg_h = [0.15, 0.05] ; %[0.08 0.08] ;
            marg_w = [0.05, 0.01] ;
        end
        
        
        % axis limits
        xlim = [min(clustNums) - 0.5, max(clustNums) + 0.5] ;
        ylim = [0, 1.05] ; 
    otherwise
        fprintf('Error: invalid plot type: %s \n', plotType)
        keyboard
end


% ------------------------------------------------------------------------
%% initialize figure window
h_main = figure('PaperPositionMode','auto','MenuBar','none',...
    'ToolBar','none','DockControls','off','Units',figUnits,...
    'OuterPosition',figPosition) ;

% --------------------------------------------------
%% loop through clusters and make plots
% NB: depends on plot type
switch plotType
    case 'pie'
        % for pie charts, need an axis for each cluster (i.e. pie chart for
        % each cluster). So loop through and make a plot
        for ind = 1:N_clusts
            % generate axis for cluster pie chart
            ax = subtightplot(n_rows, n_cols, ind, gap, ...
                marg_h, marg_w) ;
            set(ax,'Parent',h_main)
            
            % adjust axis limits
            set(ax,'xlim',xlim,'ylim',ylim)
            
            % fill graph axis
            clusterNum = ind ;
            ax = plotClusterHemilineageToAxis(ax, conn_clust_struct, ...
                clusterNum, useSubTypeFlag, titleFlag, colorMode, ...
                scaleChartFlag, sortWedgeFlag) ;
        end
        
    case 'bar'
        % for bar plots, everything goes on one axis, so make that then
        % loop through to add bars
        ax = subtightplot(1, 1, 1, gap, ...
            marg_h, marg_w) ;
        set(ax,'Parent',h_main)
        hold(ax,'on')
        
        % loop over clusters
        for ind = 1:N_clusts
            % current cluster number
            clusterNum = ind ;
            
            % add stacked bar for curent cluster
            ax = plotClusterHemilineageBarToAxis(ax, conn_clust_struct, ...
                clusterNum, useSubTypeFlag, colorMode, sortStackFlag, ...
                fullLabelFlag, horizontalFlag) ;
        end
        
        % ---------------------------------------
        % axis properties
        if horizontalFlag
            % adjust axes
            set(ax, 'xlim', ylim, 'ylim', xlim)
            set(ax, 'YTick', 1:N_clusts) ;
            
            % add x and y labels
            ylabel('interneuron cluster number')
            xlabel('hemilineage fraction')
            
            % general axis formatting
            ax = b1_paper_prettify_axis(ax, true, false, false, ...
                true, false, 0.85) ;
            
            % reverse y direction so cluster 1 is on top
            set(ax,'ydir','reverse')
        else
            % adjust axes
            set(ax, 'xlim', xlim, 'ylim', ylim)
            set(ax, 'XTick', 1:N_clusts) ;
            
            % add x and y labels
            xlabel('interneuron cluster number')
            ylabel('hemilineage fraction')
            
            % general axis formatting
            ax = b1_paper_prettify_axis(ax, true, false, false, ...
                false, true, 0.85) ;
        end
        % shrink axis ticks
        ax.TickLength = 0.3*ax.TickLength ;
    otherwise
        fprintf('Error: invalid plot type: %s \n', plotType)
        keyboard
end
% --------------------------------------------------------------
%% add key for hemilineage colors
if legendFlag && ...
        any(ismember(colorMode,{'hemilineage', 'behavior', 'neurotransmitter'}))
    % get colors for all hemilineages
    hemilineage_names = {'0A', '2A', '3B', '5B', '6A', '6B', '7B', ...
        '8B', '11A', '11B', '12A', '17A', '18B', '19A', '19B', 'abd', ...
        'emb'} ;
    [colorsCurr, colorKey] = getHemilineageColor(hemilineage_names, ...
        colorMode) ;
    
    % make sure we're only getting unique colors (if using non-hemilineage
    % cmap)
    %     [colorsCurr, ia, ~] = unique(colorsCurr,'rows','stable') ;
    %     colorKey = colorKey(ia) ;
    if ~strcmpi(colorMode,'hemilineage')
        [colorKey, ia, ~] = unique(colorKey) ;
        colorsCurr = colorsCurr(ia,:) ;
    end
    
    % with unique colors/keys, make sure "unknown" is at bottom, if it
    % exists
    unknown_idx = cellfun(@(y) strcmpi(y,'unknown'), colorKey) ;
    if sum(unknown_idx) == 1
        % to make sure it gets to the bottom, split both colors and keys
        % into the known entries (temp) and unknown
        colorsTemp = colorsCurr(~unknown_idx,:) ;
        colorUnknown = colorsCurr(unknown_idx,:) ;
        
        keysTemp = colorKey(~unknown_idx) ;
        keyUnknown = colorKey(unknown_idx) ;
        
        % ... then recombine them as top and bottom
        colorsCurr = vertcat(colorsTemp, colorUnknown) ;
        colorKey = vertcat(keysTemp, keyUnknown) ;
    end
    % create axis for putting key in
    ax_key = subtightplot(n_rows, n_cols, N_clusts+1, gap, ...
        marg_h, marg_w) ;
    
    xlim = [0,1] ;
    ylim = [0,1] ;
    
    set(ax_key, 'xlim', xlim, 'ylim', ylim)
    
    % add text objects with hemilineage color -- going to try for two
    % columns of just hemilineage name
    col_pos = linspace(xlim(1), xlim(2), 5) ;
    % switch number of key columns depending on how many colors we have
    if length(colorKey) < 7
        col_pos = col_pos(2:end-3) ;
    else
        col_pos = col_pos(2:end-2) ;
    end
    n_col_key = length(col_pos) ;
    
    row_pos = linspace(ylim(1), ylim(2), ...
        ceil(length(colorKey)./n_col_key)+4) ;
    row_pos = fliplr(row_pos(3:end-2)) ; % keep things in ascending order
    
    % add patch as background for key -- light gray to make colors a little
    % easier to see
    if length(colorKey) < 7
        bb_vert = [0.20, 0.1 ; 0.85, 0.1 ; 0.85, 0.9 ; 0.20, 0.9] ;
    else
        bb_vert = [0.20, 0.1 ; 0.75, 0.1 ; 0.75, 0.9 ; 0.20, 0.9] ;
    end
    h_p = patch(bb_vert(:,1), bb_vert(:,2), 0.6*[1,1,1], ...
        'LineStyle','none', 'FaceAlpha', 0.3) ;
    
    %     % flip color key if not using hemilineage
    %     if ~strcmpi(colorMode, 'hemilineage')
    %        colorKey = flipud(colorKey) ;
    %        colorsCurr = flipud(colorsCurr) ;
    %     end
    
    % loop through text positions
    counter = 1 ;
    for cc = 1:length(col_pos)
        for rr = 1:length(row_pos)
            % make sure we haven't run out of colors
            if counter > size(colorsCurr,1)
                continue
            end
            txt_curr = sprintf('{{\\color[rgb]{%f, %f, %f} %c}} %s', ...
                colorsCurr(counter,1), colorsCurr(counter,2), ...
                colorsCurr(counter,3), char(9632),  ...
                colorKey{counter}) ;
            text(col_pos(cc), row_pos(rr), txt_curr,...
                'FontName','arial','FontSize',6)
            counter = counter + 1 ;
        end
    end
    
    % turn off key axis
    axis(ax_key, 'off')
end
% --------------------------------------------------------------
%% save plot?
if saveFlag
    set(h_main, 'Renderer','painters')
    saveFn = ['cluster_hemilineages_' plotType] ;
    if horizontalFlag && strcmp(plotType,'bar')
        saveFn = [saveFn '_horz'] ;
    end
    if useSubTypeFlag
        saveFn = [saveFn '_withSubType'] ;
    end
    if strcmpi(colorMode,'behavior') || ...
            strcmpi(colorMode,'neurotransmitter')
        saveFn = [saveFn '_' colorMode] ;
    end
    exportgraphics(h_main, fullfile(savePath, [saveFn '.png']),...
        'Resolution', 500)
    print(h_main, fullfile(savePath, [saveFn '.svg']),'-dsvg')
end