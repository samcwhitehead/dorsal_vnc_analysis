% -------------------------------------------------------------------------
% script to take a given cluster of vnc neurons, and plot them on a
% template vnc to demonstrate what they look like
% -------------------------------------------------------------------------
%% path info and params
[mfilePath, ~, ~] = fileparts(mfilename('fullpath')) ; 
figDirectory = fileparts(mfilePath) ;
parentDirectory = fileparts(figDirectory) ; 
dataRoot = fullfile(parentDirectory, 'data') ;
saveRoot = fullfile(mfilePath, 'output', 'cluster_overlap') ; 

% data types to look at overlap for
dataName1 = 'IN' ;
dataName2 = 'DN' ;

% which DN to plot overlap for
DN_name_list = {'DNg02', 'DNpIP10', 'DNp01'} ;

% which cluster to look at
clustNumList = [10, 10, 11] ;

% how many cluster neurons to plot per DN
N_plots = 5 ;

% SAVE?
saveFlag = true ;

% how to normalize overlap data
normType = 'sym_size' ;

% use masked overlap?
maskedFlag = true ;

% params for neuron/mask data
inputDataType = 'coords_non_sym' ;
outputDataType = 'coords_sym' ;

% -------------------------------
% plot params
% fig and axis positions
fontSize = 6 ;

clusterColors = brewermap(27, 'Dark2') ;
DN_color = [0, 0, 0] ;

drawVNCFlag = true ;
drawMaskFlag = false ;
drawOverlapFlag = true ;
vncAlpha = 0.05 ;

% add text with neuron names?
labelFlag = true ;

% figure window size
figPosition = [5.4583, 1.8333, 0.52*6.2917, 0.52*7.8125] ;
% -------------------------------------
%% find associated overlap calculations
% paths to data and overlap
overlapPath = fullfile(dataRoot,'overlap_calculations') ;
if maskedFlag
    overlapPath = fullfile(overlapPath, 'masked') ;
end
dataPath1 = fullfile(dataRoot,dataName1) ;
dataPath2 = fullfile(dataRoot,dataName2) ;

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

% load the voxel count for each binarized neuron (used for normalization)
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

% --------------------------------------------------------
%% load IN clustering results
clusterPathIN = fullfile(dataRoot, 'clustering', ...
    'conn_clust_struct_IN.mat') ;
conn_clust_struct = importdata(clusterPathIN) ;

% ---------------------------------------------------------------------
%% LOOP OVER DN + CLUSTER NUM COMBINATIONS
for k = 1:length(DN_name_list)
    % get current DN, IN cluster num, savePath, and IN color
    DN_name = DN_name_list{k} ;
    clustNum = clustNumList(k) ;
    savePath = fullfile(saveRoot, ['cluster_' num2str(clustNum,'%02d') ...
        '_' DN_name ]) ;
    
    if ~exist(savePath,'dir') && saveFlag
        mkdir(savePath)
    end
    
    IN_color = clusterColors(clustNum,:) ;
    
    % ---------------------------------------------------------------------
    %% get overlap data for this cluster + DN
    % get row (IN) and col (DN) labels
    row_labels = get_neuron_names_from_overlap(overlap_data_norm, 1) ;
    col_labels = get_neuron_names_from_overlap(overlap_data_norm, 2) ;
    
    % use col labels to find DN-matching column
    col_match_idx = cellfun(@(y) strcmpi(DN_name, y), col_labels) ;
    overlap_data_curr = overlap_data_norm(:, col_match_idx) ;
    
    % sort rows according to In clusters, then pull out rows matching current
    % cluster
    [row_labels, sort_ind] = sort_neuron_labels_overlap(row_labels, ...
        'IN', true, [], [], conn_clust_struct) ;
    
    % also sort overlap data
    overlap_data_curr = overlap_data_curr(sort_ind,:) ;
    
    % finally, get cluster boundaries for shading
    [T_clust_perm, T_clust_perm_unique, ~] = ...
        get_cluster_outperm(conn_clust_struct) ;
    
    cluster_idx = (T_clust_perm == T_clust_perm_unique(clustNum)) ;
    
    overlap_data_curr = overlap_data_curr(cluster_idx,:) ;
    row_labels_curr = row_labels(cluster_idx) ;
    
    % ---------------------------------------------------------------------
    %% sort interneurons in cluster by their overlap with DN
    [overlap_data_sorted, overlap_sort_ind] = ...
        sort(overlap_data_curr{:,:},'descend') ;
    
    row_labels_sorted = row_labels_curr(overlap_sort_ind) ;
    N_plots_curr = min([N_plots, numel(row_labels_sorted)]) ;
    
    % ----------------------------------------------------------------------
    %% loop through and make neuron plots
    % initialize storage for figure windows
    fig_array = gobjects(N_plots_curr, 1) ;
    
    % always using the same DN here
    output_neuron = {DN_name, 'DN', outputDataType} ;
    
    % run loop
    for n = 1:N_plots_curr
        % print update
        fprintf('Drawing neuron ... %d/%d \n',n, N_plots_curr)
        
        % get current IN info
        input_neuron = {row_labels_sorted{n}, 'IN', inputDataType} ;
        
        % make plot
        fig_array(n) = vizMaskedOverlapFunc(input_neuron, output_neuron, ...
            IN_color, DN_color, vncAlpha, drawVNCFlag, drawMaskFlag, ...
            drawOverlapFlag, labelFlag, saveFlag, savePath, figPosition) ;
        
        % print update
        fprintf('Completed drawing neuron %d/%d \n',n, N_plots_curr)
        
        % close figure?
        if saveFlag
            close(fig_array(n))
        end
    end
end

