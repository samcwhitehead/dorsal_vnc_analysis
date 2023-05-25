% -------------------------------------------------------------------------
% FUNCTION to generate normalized overlap matrix plots
%
%{
dataName1 = 'VUM' ;
dataName2 = 'VUM' ;

h_overlap = makeNormOverlapPlotFunc(dataName1, dataName2) ; 
%}
%
% -------------------------------------------------------------------------
function [h_overlap, ax, row_ind, col_ind] = ...
    makeNormOverlapPlotFunc(dataName1, dataName2, savePath, dataRoot, ...
        saveFlag, normType, maskedFlag, threshVal, colorTickLabelsFlag,...
        clim, figPosition, matSortCase, row_sort_ind, col_sort_ind)
% ---------------------------------
%% inputs and params
if ~exist('savePath','var') || isempty(savePath)
    savePath = pwd ;
end
if ~exist('dataRoot','var') || isempty(dataRoot)
    [mfilePath, ~, ~] = fileparts(mfilename('fullpath')) ;
    figDirectory = fileparts(mfilePath) ;
    parentDirectory = fileparts(figDirectory) ;
    dataRoot = fullfile(parentDirectory, 'data') ;
end
if ~exist('saveFlag','var') || isempty(saveFlag)
    % save plot?
    saveFlag = true ;
end
if ~exist('normType','var') || isempty(normType)
    % how to normalize overlap data?
    normType = 'none' ; % 'sym_size' ; | 'none'
end
if ~exist('maskedFlag','var') || isempty(maskedFlag)
    % use masked overlap?
    maskedFlag = true ;
end
if ~exist('threshVal','var') || isempty(threshVal)
    % minimum voxel overlap count 
    % NB: quick back of the envelope says our voxels are ~0.46 micron per
    % side. if we assume a synapse is roughly 3^3 micron^3, then a
    % reasonable threshold value is around 300
    threshVal = 0 ;
end
if ~exist('colorTickLabelsFlag','var') || isempty(colorTickLabelsFlag)
    % color x and y tick labels with associated cluster label, where
    % possible?
    colorTickLabelsFlag = true ;
end
if ~exist('clim','var') || isempty(clim)
    clim = [] ;
end
if ~exist('figPosition','var') || isempty(figPosition)
    figPosition = [0.5, 0.5, 11.0, 11.0] ; 
end
if ~exist('matSortCase','var') || isempty(matSortCase)
    % how should we sort rows/cols of matrix
    matSortCase = 'new_cluster' ; % 'new_cluster' | 'io_cluster' | 'manual' | 'file_list'
end
if ~exist('row_sort_ind','var') || isempty(row_sort_ind)
    row_sort_ind = [] ; 
end
if ~exist('col_sort_ind','var') || isempty(col_sort_ind)
    col_sort_ind = [] ; 
end

% ----------------------------------
% other params
% ----------------------------------
% take table transpose?
% NB: should probably skip this with masked overlap, since now row and col
% are not interchangeable
transposeFlag = false ;

% z score by row
zScoreFlag = false ;

% ---------------------------------------------------------------
%% how to sort rows/cols of overlap matrix plot
% options for 'new_cluster' clustering 
distance_type = 'correlation' ;
if all(strcmp(dataName1,dataName2)) 
    linkage_type = 'average' ;
    clust_dim = 'row' ;
else
    linkage_type = 'average' ;  % 'single' ; %'average' ; %'complete' ;
    clust_dim = 'all' ;
end

% use the connectivity-based IN/DN clustering to organize the overlap plot?
% NB: we can hard-code this in because we have to pick the 'io_cluster'
% case for matSortCase for this to even be applied
useINClusterFlag = true ;
useDNClusterFlag = true ;

% need to load the connectivity clusters, if using
if strcmp(matSortCase, 'io_cluster')
    % IN input/output clustering
    clusterPath = fullfile(dataRoot, 'clustering', ...
        'conn_clust_struct_I*.mat') ;
    clustDir = dir(clusterPath) ; 
    if length(clustDir) ~= 1
        fprintf('Error finding IN clustering results \n') 
        keyboard
    end
    conn_clust_struct_IN = importdata(fullfile(clustDir(1).folder, ...
        clustDir(1).name)) ;
    
    % DN input/output clustering
    clusterPathDN = fullfile(dataRoot, 'connectivity', ...
        'conn_clust_struct_DN.mat') ;
    conn_clust_struct_DN = importdata(clusterPathDN) ;
    
else
    % otherwise just leave them as empty place holders
    conn_clust_struct_IN = [] ;
    conn_clust_struct_DN = [] ;
end

% load neuronColorStruct if we're supposed to color tick labels
if colorTickLabelsFlag
   colorStructPath = fullfile(dataRoot, 'connectivity', ...
       'neuronColorStruct.mat') ;
   neuronColorStruct = importdata(colorStructPath) ; 
end

% ---------------------------------------------------------
%% plot params
% cmap = magma() ;
cmap = viridis() ;
cmap = cmap(50:end,:) ; 

% -------------------------------------
%% load associated overlap calculations
% NB: with masked data, we don't have symmetry, so IN_DN ~= DN_IN
[overlap_data, vox_size_table1, vox_size_table2] = ...
    loadCombinedOverlapMat(dataName1, dataName2, dataRoot, maskedFlag) ; 


% ----------------------------------------
%% normalize overlap mat
if maskedFlag
    overlap_data_norm = normalize_masked_overlap_data(overlap_data, ...
        vox_size_table1, vox_size_table2, normType) ;
else
    % normalize overlap matrix
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

% -----------------------------------------
%% get labels
label_cell_1 = get_neuron_names_from_overlap(overlap_data_norm, 1) ;
label_cell_2 = get_neuron_names_from_overlap(overlap_data_norm, 2) ;

if strcmp(dataName1, 'MN')
    label_cell_1 = parse_mn_names(label_cell_1) ; 
end
if strcmp(dataName2, 'MN')
    label_cell_2 = parse_mn_names(label_cell_2) ; 
end
% -----------------------------------------
%% sort labels (and data)?
switch matSortCase
    case 'io_cluster'
        % in this case, use already-calculated IN and DN clustering
        % ROW labels
        if all(strcmpi(dataName1,'DN'))
            [~, row_ind] = sort_neuron_labels_overlap(label_cell_1, ...
                dataName1, useINClusterFlag, useDNClusterFlag, [], [], ...
                conn_clust_struct_DN) ;
        else
            [~, row_ind] = sort_neuron_labels_overlap(label_cell_1, ...
                dataName1, useINClusterFlag, useDNClusterFlag, [], [], ...
                conn_clust_struct_IN) ;
        end

        % COLUMN labels
        if all(strcmpi(dataName2,'DN'))
            [~, col_ind] = sort_neuron_labels_overlap(label_cell_2, ...
                dataName2, useINClusterFlag, useDNClusterFlag, [], [], ...
                conn_clust_struct_DN) ;
        else
            [~, col_ind] = sort_neuron_labels_overlap(label_cell_2, ...
                dataName2, useINClusterFlag, useDNClusterFlag, [], [], ...
                conn_clust_struct_IN) ;
        end
    
    case 'new_cluster' 
        % in this case, calculate new sorting based on hierarchical
        % clustering of overlap vals
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
        
    case 'manual'
        % in this case, use input user input variables "row_sort_ind" and
        % "col_sort_ind" to sort matrix. but first check if they're empty
        if isempty(row_sort_ind)
           row_sort_ind = 1:size(overlap_data_norm,1) ; 
        end
        if isempty(col_sort_ind)
           col_sort_ind = 1:size(overlap_data_norm,2) ;  
        end
        
        % just set the sorting indices we'll apply equal to the inputs
        row_ind = row_sort_ind ; 
        col_ind = col_sort_ind ; 
    
    case 'file_list'
        % somewhere, we should have a list of the neurons stored in an
        % order we kind of want. in this case, use that!
        useINClusterFlag = false ; 
        useDNClusterFlag = false ;
        
         [~, row_ind] = sort_neuron_labels_overlap(label_cell_1, ...
                dataName1, useINClusterFlag, useDNClusterFlag, [], [], ...
                conn_clust_struct_IN) ;
        [~, col_ind] = sort_neuron_labels_overlap(label_cell_2, ...
                dataName2, useINClusterFlag, useDNClusterFlag, [], [], ...
                conn_clust_struct_IN) ;
    otherwise
        fprintf('Invalid selection for matrix sorting \n')
        keyboard
end

% ----------------------------------------
% whichever case we used, apply sorting to overlap table and labels
label_cell_2 = label_cell_2(col_ind) ;
label_cell_1 = label_cell_1(row_ind) ;
overlap_data_norm = overlap_data_norm(row_ind, col_ind) ;


% ------------------------------------------------------
%% z score rows
if zScoreFlag
    overlap_mean = mean(overlap_data_norm{:,:},2) ;
    overlap_std = std(overlap_data_norm{:,:},0,2) ;
    overlap_data_norm{:,:} = ...
        (overlap_data_norm{:,:} - overlap_mean)./overlap_std ;
end

% ---------------------------------------------------------
%% switch rows and columns?
if transposeFlag
    % take table transpose
    overlap_data_norm = myTableTranspose(overlap_data_norm) ;
    
    % switch labels
    temp = label_cell_1 ;
    label_cell_1 = label_cell_2 ;
    label_cell_2 = temp ;
    
end

% -------------------------------------------------------------------------
%% generate plot
[h_overlap, ax] = plot_overlap_data(overlap_data_norm, label_cell_2,  ...
    label_cell_1, [], [], cmap) ;

% generate figure position if not specified 
% NB: to do this, we need to make the plot none visible
if ~isnumeric(figPosition)
    % starting position and scale of axes
    figPosition = [0.5, 0.5, 6.5, 11.0] ;
    inchPerNeuronX = 0.1 ; % 0.1
    inchPerNeuronY = 0.1 ; % 0.1
    
    %  make figure invisible
    set(h_overlap, 'visible', 'off', 'Units', 'inches') ; 
    
    % try to increase height of axis
    if strcmp(dataName1,dataName2)
        axis(ax,'equal')
    elseif ismember(dataName1, {'IN','DN'})
        ax_pos = get(ax,'Position') ;
        ax_pos(4) = ax_pos(4) + 0.065 ;
        set(ax,'Position',ax_pos)
    end
    
    % get new figure size
    n_rows = length(get(ax,'YTick')) ; 
    n_cols = length(get(ax,'XTick')) ; 
    
    figPositionNew = figPosition ; 
    figPositionNew(3) = inchPerNeuronX*n_cols ; 
    figPositionNew(4) = inchPerNeuronY*n_rows + 1.5 ; 
    
    set(h_overlap, 'OuterPosition', figPositionNew) ; %, ...
        %'PaperPosition', figPositionNew) ; 
    
    % if we're increasing figure size dramatically, try to reduce font
    if (figPositionNew(3) > 15)
       ax.XAxis.FontSize = 4.0 ;
    end
    if (figPositionNew(4) > 15)
       ax.YAxis.FontSize = 4.0 ;
    end
    
else
    % set figure size so we have matching squares
    set(h_overlap,'Units', 'inches', 'OuterPosition',figPosition)

end

%set(ax,'xlim',xlim,'ylim',ylim)

% remove axis lines
ax.XRuler.Axle.LineStyle = 'none' ;
ax.YRuler.Axle.LineStyle = 'none' ;

% set colorbar limints
if ~isempty(clim)
    caxis(clim) ;
else
    overlap_vals = overlap_data_norm{:,:} ;
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
    c.Label.String = 'volume overlap (vox num)' ;
end

% add title to clarify what we're plotting
if iscell(dataName1)
    dataName1_title = strjoin(dataName1,',') ; 
else
    dataName1_title = dataName1 ; 
end
if iscell(dataName2)
    dataName2_title = strjoin(dataName2,',') ; 
else
    dataName2_title = dataName2 ; 
end

if maskedFlag && transposeFlag
    title_str = sprintf('%s %c %s overlap', dataName1_title, ...
        char(hex2dec('2192')), dataName2_title) ; 
elseif maskedFlag && ~transposeFlag
     title_str = sprintf('%s %c %s overlap', dataName2_title, ...
         char(hex2dec('2192')), dataName1_title) ; 
else
     title_str = sprintf('%s and %s overlap', dataName2_title, ...
         dataName1_title) ; 
end
title(title_str, 'fontName','arial','fontSize', 8, 'fontWeight','normal')

% -----------------------------------------------------------------------
%% try to color axis tick labels?
if colorTickLabelsFlag
    ax = colorNeuronTickLabels(ax, neuronColorStruct) ; 
end

% ----------------------------------------
%% save results?
if saveFlag
    suffixStr = '' ;
    if zScoreFlag
        suffixStr = [suffixStr '_zScore'] ;
    end
    if strcmpi(matSortCase, 'new_cluster')
        suffixStr = [suffixStr '_sort'] ;
    end
    %set(h_overlap,'OuterPosition',figPosition)
    savefig(h_overlap, fullfile(savePath, [dataName1 '_' dataName2 ...
        '_' normType '_overlap' suffixStr '.fig'])) ;
    exportgraphics(h_overlap, fullfile(savePath, [dataName1 '_' dataName2 ...
        '_' normType  '_overlap' suffixStr '.png']),'Resolution',600) ;
    exportgraphics(h_overlap, fullfile(savePath, [dataName1 '_' dataName2 ...
        '_' normType  '_overlap' suffixStr '.jpg']),'Resolution',600) ;
%     export_fig(h_overlap, fullfile(savePath, [dataName1 '_' dataName2 ...
%         '_' normType  '_overlap' suffixStr '.png']),'-dpng','-r600') ;
    %     print(h_overlap, fullfile(savePath, [dataName1 '_' dataName2 ...
    %         '_' normType  '_overlap' suffixStr '.svg']),'-dsvg') ;
end
end