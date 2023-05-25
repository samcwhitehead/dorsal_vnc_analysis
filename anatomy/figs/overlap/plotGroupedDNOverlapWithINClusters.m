% -------------------------------------------------------------------------
% quick and dirty script to make plots that show IN overlap with DNs. INs
% (on y axis) GROUPED according to i/o clustering
% -------------------------------------------------------------------------
% -----------------------
%% path info and params
[mfilePath, ~, ~] = fileparts(mfilename('fullpath')) ; 
figDirectory = fileparts(mfilePath) ;
parentDirectory = fileparts(figDirectory) ; 
dataRoot = fullfile(parentDirectory, 'data') ;

% data types to look at overlap for
dataName1 = 'IN' ;
dataName2 = 'DN' ;

% colormap for heatmap
cmap = viridis() ;
cmap = cmap(50:end,:) ;

% SAVE?
% saveFlag1 = false ; % bar plots
saveFlag = false ; % heatmap
savePath = fullfile(mfilePath, 'output', 'cluster_overlap') ; 
if ~exist(savePath,'dir') && saveFlag
    mkdir(savePath)
end

% how to normalize overlap data
normType = 'none' ;

% use masked overlap?
maskedFlag = true ;

% threshold overlap values?
threshVal = 300 ; % units of voxel count

% sort matrix rows and columns?
sortMatFlag = true ; 

% color y tick labels (INs)?
colorTickLabelsFlag = false ;

% -------------------------------
% plot params
% fig and axis positions
figPosition = [0.5, 2.0, 6.5, 4.0] ;
figUnits = 'inches' ;

gap = [0.00, 0.00] ; %[0.1 0.10] ;
% marg_h = [0.05, 0.00] ; % [0.05, 0.01] %[0.1 0.02] ;
marg_h = [0.1, 0.05] ;
marg_w = [0.22, 0.10] ;

% line properties
lineWidthThick = 1.0 ;
lineWidthThin = 0.5 ;
lineAlpha = 0.25 ;

% axis font size to keep neuron labels from getting crammed
fontSizeSmall = 4.25 ;

% axis limits
xlim = [] ;  % leaving it empty lets it go auto
clim = [] ;
% ylim takes care of itself
% -------------------------------------
%% find associated overlap calculations
% paths to data and overlap
overlapPath = fullfile(dataRoot,'overlap_calculations') ;
if maskedFlag
    overlapPath = fullfile(overlapPath, 'masked') ;
end

% try different permutations of the two names
overlapDataFilename1 = fullfile(overlapPath,[dataName1 '_' dataName2 ...
    '_overlap_table.mat']) ;
overlapDataFilename2 = fullfile(overlapPath,[dataName2 '_' dataName1 ...
    '_overlap_table.mat']) ;
if exist(overlapDataFilename1,'file')
    overlapDataFilename = overlapDataFilename1 ;
elseif ~exist(overlapDataFilename1,'file') && ...
        exist(overlapDataFilename2,'file')
    overlapDataFilename = overlapDataFilename2 ;
else
    fprintf('could not locate either : %s or %s \n',overlapDataFilename1, ...
        overlapDataFilename2)
end

% -------------------------------------------------------------------------
%% load overlap data
overlap_data = importdata(overlapDataFilename) ;

vox_size_table1 = importdata(fullfile(dataRoot,'vox_size_tables', ...
    dataName1, 'vox_size_table_masked.mat')) ;
vox_size_table2 = importdata(fullfile(dataRoot,'vox_size_tables', ...
    dataName2, 'vox_size_table_masked.mat')) ;

% number of neurons according to IN directory
N_stacks1 = size(vox_size_table1,1) ;
N_stacks2 = size(vox_size_table2,1) ;

% see if our matrix dimensions are off
if (size(overlap_data,1) == N_stacks2) && ...
        (size(overlap_data,2) == N_stacks1) && (N_stacks1 ~= N_stacks2)
    % if so, take transpose of overlap data
    overlap_data = myTableTranspose(overlap_data) ;
    
end

% ----------------------------------------
%% normalize overlap mat
if maskedFlag
    overlap_data_norm = normalize_masked_overlap_data(overlap_data, ...
        vox_size_table1, vox_size_table2, normType) ;
else
    overlap_data_norm = normalize_overlap_data(overlap_data, ...
        vox_size_table1, vox_size_table2, normType) ;
end

% -----------------------------------------
%% threshold overlap values?
% get indices where overlap voxel count is below threshold
overlap_vals_thresh = overlap_data_norm{:,:} ;
below_thresh_idx = (overlap_vals_thresh <= threshVal) ;


% set below thresh entries to zero
overlap_vals_thresh(below_thresh_idx) = 0 ;

% put thresholded values back in table
overlap_data_norm{:,:} = overlap_vals_thresh ;

% ---------------------------------------------
%% get row and column labels
row_labels = get_neuron_names_from_overlap(overlap_data_norm, 1) ;
col_labels = get_neuron_names_from_overlap(overlap_data_norm, 2) ;

% --------------------------------------------------------
%% load IN clustering results
clusterPathIN = fullfile(dataRoot, 'clustering', ...
    'conn_clust_struct_IN.mat') ;
conn_clust_struct_IN = importdata(clusterPathIN) ;

% ------------------------------------------
%% sort rows according to IN clusters
[row_labels, row_sort_ind] = sort_neuron_labels_overlap(row_labels, ...
    'IN', true, [], [], conn_clust_struct_IN) ;

overlap_data_norm = overlap_data_norm(row_sort_ind, :) ;

% -------------------------------------------------------------
%% group overlap matrix data by IN cluster number
% get cluster boundaries to loop over
[T_clust_perm, T_clust_perm_unique, ~] = ...
    get_cluster_outperm(conn_clust_struct_IN) ;
N_clust = length(T_clust_perm_unique) ;

% initialize new table to store output
varTypes = cell(size(col_labels)) ;
varTypes(:) = {'double'} ;

overlap_grouped = table('Size', [N_clust, size(overlap_data_norm,2)],...
    'VariableTypes', varTypes, 'VariableNames', col_labels) ;
row_labels_grouped = cell(N_clust,1) ;

exclude_idx = false(N_clust,1) ; 
for n = 1:N_clust
    % get number & index for current cluster
    clustCurr = T_clust_perm_unique(n) ;
    cluster_idx = (T_clust_perm == clustCurr) ;
    if sum(cluster_idx) <= 2
       exclude_idx(n) = true ;  
    end
        
    % take median over this range and add to grouped table
    overlap_grouped{n,:} = nanmedian(overlap_data_norm{cluster_idx,:},1) ;
%     overlap_grouped{n,:} = nanmean(overlap_data_norm{cluster_idx,:},1) ;
    
    % just store row number as string in cell array for convenience
    row_labels_grouped{n} = num2str(n) ;
end

% --------------------------------------------------------------
%% remove clusters with small neuron number 
overlap_grouped = overlap_grouped(~exclude_idx, :) ; 
row_labels_grouped = row_labels_grouped(~exclude_idx); 

% --------------------------------------------------------------
%% sort matrix rows and columns to highlight patterns?
if sortMatFlag
    % options for sorting
    distance_type = 'correlation' ;
    if all(strcmp(dataName1,dataName2))
        linkage_type = 'average' ;
        clust_dim = 'row' ;
    else
        linkage_type = 'average' ;  % 'single' ; %'average' ; %'complete' ;
        clust_dim = 'all' ;
    end
    
    % get sort ind
    [~, row_ind, col_ind] = myMatrixClusterSort(overlap_grouped{:,:},...
        clust_dim, distance_type, linkage_type) ;
    
    % apply sort ind
    overlap_grouped = overlap_grouped(row_ind, col_ind) ; 
    row_labels_grouped = row_labels_grouped(row_ind) ; 
    col_labels = col_labels(col_ind) ; 
    
end
% --------------------------------------------------------------
%% plot grouped overlap
% initialize figure window
[h_overlap, ax] = plot_overlap_data(overlap_grouped, col_labels,  ...
    row_labels_grouped, [], [], cmap) ;

%  [h_overlap, ax] = plot_overlap_data(overlap_data, labels_x, ...
%     labels_y, plotStyle, normFlag, cmap)
% set figure size so we have matching squares
set(h_overlap,'Units', 'inches', 'OuterPosition',figPosition)

%set(ax,'xlim',xlim,'ylim',ylim)

% remove axis lines
ax.XRuler.Axle.LineStyle = 'none' ;
ax.YRuler.Axle.LineStyle = 'none' ;

% add y label
ylabel(ax, 'interneuron cluster', 'FontName','arial', 'FontSize', 6)

% set colorbar limints
if ~isempty(clim)
    caxis(clim) ;
else
    overlap_vals = overlap_grouped{:,:} ;
    cmax = prctile(overlap_vals(:),99) ;
    if cmax == 0
        cmax = max(overlap_vals(:)) ;
    end
    caxis([0, cmax]) ;
end

% adjust colorbar label
if strcmpi(normType,'none')
    %     test = findobj(ax) ;
    c = findobj(h_overlap, 'Type', 'colorbar') ;
    c.Label.String = 'median volume overlap (vox num)' ;
end

% add title to clarify
if maskedFlag
    title_str = sprintf('%s %c clustered %s overlap', dataName2, ...
        char(hex2dec('2192')), dataName1) ;
else
    title_str = sprintf('%s and %s overlap', dataName2, dataName1) ;
end
title(title_str, 'fontName','arial','fontSize', 8, 'fontWeight','normal')

% ------------------------------------------------------------
%% color tick labels (just DNs here)
if colorTickLabelsFlag
    % load neuron color struct
    colorStructPath = fullfile(dataRoot, 'neuronColorStruct.mat') ;
    neuronColorStruct = importdata(colorStructPath) ;
    
    % make sure axis tick label interpreter is 'Tex'
    ax.TickLabelInterpreter = 'tex' ;
    
    % start with x ticks. loop over labels. if neuron label is in our
    % struct, change the tick label to that color
    x_tick_labels = get(ax, 'XTickLabel') ;
    for xn = 1:length(x_tick_labels)
        % current tick label
        tick_label_curr = strrep(x_tick_labels{xn},'-','_') ;
        
        % check that we have a color match
        if isfield(neuronColorStruct, tick_label_curr)
            % if we find a matching color, read it out and replace tick label
            neuronColor = neuronColorStruct.(tick_label_curr) ;
            tick_label_new = sprintf('{{\\color[rgb]{%f, %f, %f} %s}}', ...
                neuronColor(1), neuronColor(2), neuronColor(3), ...
                tick_label_curr) ;
        else
            % if we don't find a match, just get it to print black?
            tick_label_new = sprintf('{{\\color[rgb]{%f, %f, %f} %s}}', ...
                0.0, 0.0, 0.0, tick_label_curr) ;
        end
        
        % add label back to array; make sure to switch out hyphen and
        % underscore
        x_tick_labels{xn} = strrep(tick_label_new,'_','-') ;
        
    end
    
    % set axis to have new x tick labels
    set(ax,'XTickLabel',x_tick_labels) ;
end


% --------------------------------------------------------
%% save plot?
if saveFlag
    saveNameFull = fullfile(savePath, ['grouped_IN_overlap_heatmap']) ;
    savefig(h_overlap, [saveNameFull '.fig'])
    export_fig(h_overlap, saveNameFull, '-dpng','-r600')
    print(h_overlap, saveNameFull, '-dsvg')
end


% % color bar position
% cpos = c.Position ;
% cpos(4) = cpos(4)/3 ; % shrink height
% % cpos(2) = cpos(2) + 0.25 ; % shift up?
% cpos(1) = cpos(1) + 0.2 ; % shift right?
% c.Position = cpos ;
%
% % color bar label position
% clabel_pos = c.Label.Position ;
% clabel_pos(1) = clabel_pos(1) - 3.5 ;
% c.Label.Position = clabel_pos ;
