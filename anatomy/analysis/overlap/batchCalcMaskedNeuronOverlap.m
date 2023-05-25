% -------------------------------------------------------------------------
% function to perform batch calculation of masked overlap between sets of
% neurons (i.e. compare all INs against all DNs)
%
% with the introduction of input/output neuropil masks, we now have a loss
% of symmmetry -- there is a difference between looking at the overlap
% between (inputs of neuron 1 vs outputs of neuron 2) and (outputs of
% neuron 1 vs inputs of neuron 2).
%
% in this code, we'll use the convention that an entry in the overlap table
% (just switching to tables full time now -- no more of that switching data
% structure shit) will have the form:
%   overlap_table(i,j) = overlap of inputs from neuron i with outputs of
%       neuron j
%
% may have to adjust the way that normalization takes place, but can cross
% that bridge later
%
%{

EXAMPLE USAGE:

dataRoot = 'D:\Fly Imaging\Erica Interneuron Stacks\' ;
dataType1 = 'MN' ; 
dataType2 = 'MN' ; 
savePath = fullfile(dataRoot, 'overlap_calculations','masked') ;
dataMode = 'coords' ;

overlap_data = batchCalcMaskedNeuronOverlap(dataRoot, ...
    dataType1, dataType2, savePath, dataMode) ; 
%}
% -------------------------------------------------------------------------
function overlap_data = batchCalcMaskedNeuronOverlap(dataRoot, ...
    dataType1, dataType2, savePath, dataMode)
% --------------------------------------------------------------
%% set params and get data
% where to save overlap data
if ~exist('savePath','var') || isempty(savePath)
    savePath = fullfile(dataRoot, 'overlap_calculations','masked') ;
end
% use bw images or coordinate lists for overlap calculations? coordinates
% should be faster
if ~exist('dataMode','var') || isempty(dataMode)
    dataMode = 'coords' ; %'coords' | 'bw'
end

% flag options
saveFlag = true ;
overWriteFlag = false ;
plotFlag = false ;

% misc params
INPUT = 1 ; 
OUTPUT = 2 ; 

% haven't written code to work with images yet, so don't allow this to
% continue if that's the selected option
if strcmpi(dataMode, 'bw')
    fprintf('Under construction! \n')
    keyboard
end

% also check if we're comparing same neuron type against itself 
if strcmp(dataType1, dataType2)
    sameTypeFlag = true ; 
else
    sameTypeFlag = false ; 
end
% -------------------------------------------------------------------------
%% get directories for neurons
% pretty much have to use symmetrized data, since i don't have unilateral
% DN images
fileSearchStr = ['*_' dataMode '.mat'] ;

% get directories for image data
dataDir1 = dir(fullfile(dataRoot,dataType1,'binarized_sym',fileSearchStr));
dataDir2 = dir(fullfile(dataRoot,dataType2,'binarized_sym',fileSearchStr));

% how many data files in each directory?
N_stacks1 = length(dataDir1) ;
N_stacks2 = length(dataDir2) ;

% get directories for masks
if strcmp(dataMode, 'coords')
    maskSearchStr = '*_coords.mat' ;
else
    maskSearchStr = '*_mask.mat' ;
end

maskDir1 = dir(fullfile(dataRoot,dataType1,'io_masks',maskSearchStr));
maskDir2 = dir(fullfile(dataRoot,dataType2,'io_masks',maskSearchStr));

% make sure we have as many masks as stacks/coords data files
if (length(maskDir1) ~= N_stacks1) || (length(maskDir2) ~= N_stacks2)
   fprintf('Error: mismatch between masks and image data files \n')
   keyboard
end

% also define savepath for overlap data. NB: since we're looking for
% directional connections now, we lose the symmetry in filenames, i.e.
% "IN_DN" is not equivalent to "DN_IN"
overlapFilename = fullfile(savePath,[dataType1 '_' dataType2 ...
    '_overlap_table.mat']) ;

% -------------------------------------------------------------------------
%% get neuron labels
label_cell_1 = cell(N_stacks1,1) ;
for k = 1:N_stacks1
    fn = dataDir1(k).name ;
    label_cell_1{k} = getLineName(fn, dataType1, true) ;
end

label_cell_2 = cell(N_stacks2,1) ;
for k = 1:N_stacks2
    fn = dataDir2(k).name ;
    label_cell_2{k} = getLineName(fn, dataType2, true) ;
end

% -----------------------------------------------------------------
%% check if we already made some progress on overlap calculation
if (~overWriteFlag) && exist(overlapFilename,'file')
    overlap_table = importdata(overlapFilename) ;
    
    % -----------------------------------------------------------
    %% first check if we're adding any new neurons.
    % if so, append rows/cols for them. start with rows
    curr_rows = overlap_table.Properties.RowNames ;
    new_rows = setdiff(label_cell_1, curr_rows) ;
    
    % append new rows
    if ~isempty(new_rows)
        % construct new rows into their own table, filled with nans
        new_row_vals = nan(length(new_rows), size(overlap_table,2)) ;
        new_row_table = array2table(new_row_vals) ;
        new_row_table.Properties.RowNames = new_rows ;
        new_row_table.Properties.VariableNames = ...
            overlap_table.Properties.VariableNames ;
        
        % append new row table
        overlap_table = [overlap_table ; new_row_table] ;
        
        % sort row labels/dir by this new order
        rows_all = overlap_table.Properties.RowNames ;
        [~, sort_idx1] = ismember(rows_all , label_cell_1) ;
        label_cell_1 = label_cell_1(sort_idx1) ;
        bwDir1 = dataDir1(sort_idx1) ;
    end
    
    % then find new columns
    curr_cols = overlap_table.Properties.VariableNames ;
    new_cols = setdiff(label_cell_2, curr_cols) ;
    
    % ... and append them as well
    if ~isempty(new_cols)
        % construct new rows into their own table, filled with nans
        new_col_vals = nan(size(overlap_table,1), length(new_cols)) ;
        new_col_table = array2table(new_col_vals) ;
        new_col_table.Properties.RowNames = overlap_table.Properties.RowNames ;
        new_col_table.Properties.VariableNames = new_cols ;
        
        % append new row table
        overlap_table = [overlap_table , new_col_table] ;
        
        % sort col labels/dir by this new order
        vars_all = overlap_table.Properties.VariableNames ;
        [~, sort_idx2] = ismember(vars_all , label_cell_2) ;
        label_cell_2 = label_cell_2(sort_idx2) ;
        bwDir2 = dataDir2(sort_idx2) ;
    end
    
    % -------------------------------------------
    %% then try to find where we last let off
    overlap_vals = overlap_table{:,:} ;
    idx1_start = find((overlap_vals(:) ~= 0) & ...
        (overlap_vals(:) ~= 1),1,'last') ;
    [ind1_start, ~] = ind2sub(size(overlap_vals),idx1_start) ;
else
    % -------------------------------------------------
    %% if no data exists thusfar, start new table
    temp = nan(N_stacks1, N_stacks2) ;
    overlap_table = array2table(temp) ;
    overlap_table.Properties.RowNames = label_cell_1 ;
    overlap_table.Properties.VariableNames = label_cell_2 ;
    ind1_start = 1 ;
end

% -----------------------------------------------------------
%% loop over images to fill in table
for ind1 = ind1_start:N_stacks1
    % -------------------------------------
    % data for current neuron of dataType1
    dataFilename_1 = fullfile(dataDir1(ind1).folder, dataDir1(ind1).name) ;
    lineName1 = label_cell_1{ind1} ;
    data1 = importdata(dataFilename_1) ;
    
    % masks for current neuron (type 1) -- counting on mask directory being
    % in the same order here
    maskFilename_1 = fullfile(maskDir1(ind1).folder, maskDir1(ind1).name);
    mask_cell_1 = importdata(maskFilename_1) ; 
    
    % according to our convention, we only need the INPUT mask for neurons
    % of dataType1
    mask1 = mask_cell_1{INPUT} ; 
    
    % apply mask to data1 
    data1_masked = intersect(data1, mask1) ;
    
    % make second loop over images
    for ind2 = 1:N_stacks2
        % if we're comparing a set of neurons against itself, skip diagonal
        % entries
        if sameTypeFlag && (ind1 == ind2)
            continue
        end
        
        % ... otherwise proceed with loading data for neuron of type 2
        tic
        % ----------------
        % filename and line name for current neuron of dataType2
        dataFilename_2 = fullfile(dataDir2(ind2).folder, dataDir2(ind2).name) ;
        lineName2 = label_cell_2{ind2} ;
        
        % check to see if this entry is already filled
        if ~isnan(overlap_table{lineName1, lineName2}) && ...
                ~overWriteFlag
            fprintf('Already calculated overlap for %s vs %s \n',...
                lineName1, lineName2)
            continue
        end
        
        % otherwise load data and mask for current neuron of dataType2
        data2 = importdata(dataFilename_2) ;
        
        maskFilename_2 = fullfile(maskDir2(ind2).folder, maskDir2(ind2).name);
        mask_cell_2 = importdata(maskFilename_2) ; 
        mask2 = mask_cell_2{OUTPUT} ; 
        
        data2_masked = intersect(data2, mask2) ; 
        % --------------------
        % calculate overlap
        overlap_table{lineName1, lineName2} = ...
            numel(intersect(data1_masked, data2_masked)) ;
        toc
        
    end
    % periodically save so we don't lose data if restart happens
    if saveFlag
        save(overlapFilename, 'overlap_table')
    end
end

% ----------------------------------------------------------------
%% sort table entries alphabetically
[rowNamesSort, row_sort_idx] = ...
    sort(overlap_table.Properties.RowNames) ;
[colNamesSort, col_sort_idx] = ...
    sort(overlap_table.Properties.VariableNames) ;
overlap_table = overlap_table(row_sort_idx, col_sort_idx) ;

% % ----------------------------------------------------------------
% %% fill nan entries that shouldn't be there
% % (these would arise if we append a row/column for a new neuron)
% [nan_row, nan_col] = find(isnan(overlap_table{:,:})) ;
% 
% % there are lots of intended nan entries for square overlap
% % matrices (when we compare a set against itself), so make sure to
% % only grab upper triangular elements in those cases
% if strcmp(dataName1,dataName2)
%     to_fill_idx = (nan_row < nan_col) ;
%     to_fill_row = nan_row(to_fill_idx) ;
%     to_fill_col = nan_col(to_fill_idx) ;
% else
%     to_fill_row = nan_row ;
%     to_fill_col = nan_col ;
% end
% 
% % loop through indices that need filling
% for ind = 1:length(to_fill_row)
%     tic
%     % get appropriate data file for current ROW
%     rowNameCurr = rowNamesSort{to_fill_row(ind)} ;
%     fn_idx_row = cellfun(@(y) strcmp(y, rowNameCurr), label_cell_1) ;
%     dataFilename_row = fullfile(bwDir1(fn_idx_row).folder, ...
%         bwDir1(fn_idx_row).name) ;
%     
%     % get appropriate data file for current COLUMN
%     colNameCurr = colNamesSort{to_fill_col(ind)} ;
%     fn_idx_col = cellfun(@(y) strcmp(y, colNameCurr), ...
%         label_cell_2) ;
%     dataFilename_col = fullfile(bwDir2(fn_idx_col).folder, ...
%         bwDir2(fn_idx_col).name) ;
%     
%     % load both data sets
%     bwMatRow = importdata(dataFilename_row) ;
%     bwMatCol = importdata(dataFilename_col) ;
%     
%     % get overlap and add to table
%     overlap_table{rowNameCurr, colNameCurr} = ...
%         sum(bwMatRow(:) & bwMatCol(:)) ;
%     toc
%     
% end
% %         nan_row_ind = find(all_nan_row_idx) ;
% %         for indd = nan_row_ind
% %             fn_idx = cellfun(@(y) strcmp(
% %             dataFilename1 =
% %         end

% -------------------------------------------------------------
%% save final results?
if saveFlag
    save(overlapFilename, 'overlap_table')
end

% give generic name for output
overlap_data = overlap_table ;

% ------------------------------------------------------------------------
%% plot results?
if plotFlag
    %save(fullfile(dataRoot,'overlap_mat_backup.mat'), 'overlap_mat')
    % get labels
    h_overlap = plot_overlap_data(overlap_data, label_cell_2,  label_cell_1) ;
    
    if saveFlag
        % make figure full-screen
        set(h_overlap, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
        print(h_overlap, fullfile(savePath, [dataName1 '_' dataName2 ...
            '_overlap.png']), '-dpng','-r300')
    end
    
end
end