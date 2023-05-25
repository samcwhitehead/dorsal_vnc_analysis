% -------------------------------------------------------------------------
% script to make a tree plot showing the INs that a given DN connects to 
% -------------------------------------------------------------------------
%% path and param info
[mfilePath, ~, ~] = fileparts(mfilename('fullpath')) ; 
figDirectory = fileparts(mfilePath) ;
parentDirectory = fileparts(figDirectory) ; 
dataRoot = fullfile(parentDirectory, 'data') ;
overlapPath = fullfile(dataRoot, 'overlap_calculations') ;

savePath = fullfile(mfilePath, 'output', 'tree_plots') ;
if ~exist(savePath, 'dir')
    mkdir(savePath)
end

% parent neuron name
dataName1 = 'IN' ; 
dataName2 = 'DN' ; 
parentNode = 'DNpIP10' ; % 'DNpIP10' | 'DNp01' | 'DNp07' | 'DNMW' | 'DNg02' | ...

% data params
saveFlag = true ; % save plots
normType = 'none' ; % how to normalize overlap data? % 'sym_size' ; | 'none'
maskedFlag = true ; % use masked overlap?
threshVal = 300 ; % minimum value of overlap. 500 ~ a 3.7x3.7x3.7 micron cube
sortMatFlag = true ; % perform hierarchical clustering on overlap to get new sorting

% plot params
maxNumNodes = 20 ; 
sortEdgeFlag = true ; % sort child nodes by edge strength?
figPosition = [6.8646, 4.7917, 2.4, 2.5] ; 

cmap = viridis() ;
cmap = [[0,0,0]; cmap(50:end,:)] ;
clim = [0, 1.65e4] ; 
cmap_struct.cmap = cmap ; 
cmap_struct.clim = clim ; 

arrowColor = [0.7, 0, 0] ; 
nodeColor = [0, 0, 0] ; 

xlim = [-0.05, 1.1] ; 
ylim = [-0.05, 1.1] ; 

% make sure save path exists
if ~exist(savePath,'dir')
    mkdir(savePath)
end

% -------------------------------------------------------
%% load/process overlap data
% get overlap data table
[overlap_data, vox_size_table1, vox_size_table2] = ...
    loadCombinedOverlapMat(dataName1, dataName2, dataRoot, maskedFlag) ; 

% --------------------
% normalize 
if maskedFlag
    overlap_data_norm = normalize_masked_overlap_data(overlap_data, ...
        vox_size_table1, vox_size_table2, normType) ;
else
    % normalize overlap matrix
    overlap_data_norm = normalize_overlap_data(overlap_data, ...
        vox_size_table1, vox_size_table2, normType) ;
end

% -----------------------------------------
% threshold overlap values
% get indices where overlap voxel count is below threshold
overlap_vals_thresh = overlap_data_norm{:,:} ;
below_thresh_idx = (overlap_vals_thresh <= threshVal) ; 

% set below thresh entries to zero
overlap_vals_thresh(below_thresh_idx) = 0 ; 

% put thresholded values back in table
overlap_data_norm{:,:} = overlap_vals_thresh ; 

% -----------------------------------------
% get labels
label_cell_1 = get_neuron_names_from_overlap(overlap_data_norm, 1) ;
label_cell_2 = get_neuron_names_from_overlap(overlap_data_norm, 2) ;

% -----------------------------------------
% sort labels (and data)?
if sortMatFlag
    % in this case, calculate new sorting based on hierarchical
    % clustering of overlap vals
    % set clustering params
    distance_type = 'correlation' ;
    linkage_type = 'average' ;  % 'single' ; %'average' ; %'complete' ;
    clust_dim = 'all' ;

    % read out overlap values
    overlap_vals = overlap_data_norm{:,:} ;
    
    % sort both rows and columns (should make some of these options
    % variables...)
    [~, row_ind, col_ind] = myMatrixClusterSort(overlap_vals,...
        clust_dim, distance_type, linkage_type) ;
    
    % if square matrix, try to make the rows and cols match
    if strcmp(clust_dim,'row') && ...
            (size(overlap_data_norm,1) == size(overlap_data_norm,2))
        col_ind = row_ind ;
    end

    %  apply sorting to overlap table and labels
    label_cell_2 = label_cell_2(col_ind) ;
    label_cell_1 = label_cell_1(row_ind) ;
    overlap_data_norm = overlap_data_norm(row_ind, col_ind) ;   
end

% ----------------------------------------------------------------
%% pull out column corresponding to parent node
parent_idx = cellfun(@(y) strcmp(y, parentNode), label_cell_2) ; 
if sum(parent_idx) ~= 1
   fprintf('Trouble finding column for parent node %s ...\n',parentNode) 
end
connections = overlap_data_norm{:,parent_idx}' ; 
childNodes = label_cell_1 ; 

% ----------------------------------------------------------------
%% make tree plot
% initialize figure window
h_main =  figure('PaperPositionMode','auto', ...
    'Units', 'inches', ...
    'OuterPosition', figPosition) ;

% intialize axis
gap = [0.00, 0.00] ; %[0.1 0.10] ;
marg_h = [0.35, 0.1] ; %[0.1 0.02] ;
marg_w =  [0.00, 0.00] ; % [0.125, 0.125] ;

ax = subtightplot(1,1,1, gap, marg_h, marg_w) ;
hold on

% draw tree to axis
[ax, hg] = myPlotTreeToAx(ax, parentNode, childNodes, connections, ...
    maxNumNodes, sortEdgeFlag, cmap_struct,  nodeColor, arrowColor) ; 

% set axis limits
set(ax, 'xlim', xlim, 'ylim', ylim)

% -------------------------------------------------
%% save tree plot
if saveFlag
    % get filename to save to 
    suffixStr = '' ; 
    if sortMatFlag
        suffixStr = [suffixStr '_sort'] ; 
    end
   saveName = sprintf('%s_tree%s',parentNode, suffixStr) ; 
   
   % save png, svg, ?
   print(h_main, fullfile(savePath, [saveName '.svg']), '-dsvg') ;  
   exportgraphics(h_main, fullfile(savePath, [saveName '.png']),...
       'Resolution',600)
end