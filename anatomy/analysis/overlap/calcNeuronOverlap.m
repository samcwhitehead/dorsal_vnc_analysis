% -------------------------------------------------------------------------
% function to look at overlap of binarized confocal images of neurons.
% meant to replace the "test_neuron_overlap_v*.m" script
%{
dataRoot1 = 'D:\Fly Imaging\Erica Interneuron Stacks\MN\' ; 
dataRoot2 = 'D:\Fly Imaging\Erica Interneuron Stacks\MN\' ;
savePath = 'D:\Fly Imaging\Erica Interneuron Stacks\overlap_calculations\' ; 
overlap_data = calcNeuronOverlap(dataRoot1, dataRoot2, savePath) ;
%}
% -------------------------------------------------------------------------
function overlap_data = calcNeuronOverlap(dataRoot1, dataRoot2, savePath, ...
    dataMode) 
% --------------------------------------------------------------
%% set params and get data 
% (assumes stack categories are in different folders)
if ~exist('savePath','var') || isempty(savePath) 
   savePath = 'D:\Fly Imaging\Erica Interneuron Stacks\overlap_calculations\' ; 
end
if ~exist('dataMode','var') || isempty(dataMode) 
   dataMode = 'table' ; %'struct' | 'mat' | 'table'
end
bwDir1 = dir(fullfile(dataRoot1,'binarized_new', '*bw.mat')) ;
N_stacks1 = length(bwDir1) ;
bwDir2 = dir(fullfile(dataRoot2,'binarized_new', '*bw.mat')) ;
N_stacks2 = length(bwDir2) ;

% imHeight = 1119 ; 
% imWidth = 573 ; 
% imDepth = 219 ; 

saveFlag = true ;
overWriteFlag = false ;
plotFlag = false ;

dataRoot1_split = strsplit(dataRoot1,'\') ; 
dataName1 = dataRoot1_split{end-1} ; 
dataRoot2_split = strsplit(dataRoot2,'\') ; 
dataName2 = dataRoot2_split{end-1} ; 

overlapFilename = fullfile(savePath,[dataName1 '_' dataName2 ...
    '_overlap_' dataMode '.mat']) ;
% -------------------------------------------------------------------------
%% get neuron labels
label_cell_1 = cell(N_stacks1,1) ;
for k = 1:N_stacks1
    fn = bwDir1(k).name ;
    label_cell_1{k} = getLineName(fn, dataName1, true) ;
end

label_cell_2 = cell(N_stacks2,1) ;
for k = 1:N_stacks2
    fn = bwDir2(k).name ;
    label_cell_2{k} = getLineName(fn, dataName2, true) ;
end

% -------------------------------------------------------------------------
%% switch approach depending on how we're saving data
switch dataMode
    case 'mat'
        % -----------------------------------------------------------------
        % check if we already made some progress on overlap calculation
        if (~overWriteFlag) && exist(overlapFilename,'file')
            overlap_mat = importdata(overlapFilename) ;
            idx1_start = find((overlap_mat(:) ~= 0) & ...
                (overlap_mat(:) ~= 1),1,'last') ;
            [ind1_start, ~] = ind2sub(size(overlap_mat),idx1_start) ;
        else
            overlap_mat = nan(N_stacks1, N_stacks2) ;
            ind1_start = 1 ;
        end
        % -----------------------------------------------------------
        % loop over images
        for ind1 = ind1_start:N_stacks1
            % ---------------------
            % bw image 1
            dataFilename_1 = fullfile(bwDir1(ind1).folder, bwDir1(ind1).name) ;
            bwMat1 = importdata(dataFilename_1) ;
            
            if strcmp(dataName1,dataName2)
                ind2_start = (ind1+1) ;
            else
                ind2_start = 1 ;
            end
            % make second loop over images
            for ind2 = ind2_start:N_stacks2
                tic
                % ----------------
                % image 2
                dataFilename_2 = fullfile(bwDir2(ind2).folder, bwDir2(ind2).name) ;
                bwMat2 = importdata(dataFilename_2) ;
                
                % --------------------
                % get overlap
                overlap_mat(ind1, ind2) = sum(bwMat1(:) & bwMat2(:)) ;
                toc
                
            end
            % periodically save so we don't lose data if restart happens
            if saveFlag
                save(overlapFilename, 'overlap_mat')
            end
        end
        
        % try to save overlap matrix as table also
        try
            overlapTableFilename = fullfile(savePath,[dataName1 '_' ...
                dataName2 '_overlap_table.mat']) ;
            overlap_table = array2table(overlap_mat, 'RowNames', ...
                label_cell_1, 'VariableNames', label_cell_2) ; 
            save(overlapTableFilename, 'overlap_table')
        catch
           fprintf('Could not save to table \n') 
        end
        overlap_data = overlap_mat ; 
        
    case 'struct'
        % -------------------------------------------------------------------------
        % check if we already made some progress on overlap calculation
        if (~overWriteFlag) && exist(overlapFilename,'file')
            overlap_struct = importdata(overlapFilename) ;
        else
            overlap_struct = struct('label_1',[],'label_2',[],...
                'overlap_val',[]) ;
        end
        % -----------------------------------------------------------
        % loop over images
        for ind1 = 1:N_stacks1
            % ---------------------
            % get filename for first image
            dataFilename_1 = fullfile(bwDir1(ind1).folder, bwDir1(ind1).name) ;
            
            % determine starting index for inner loop
            if strcmp(dataName1,dataName2)
                ind2_start = (ind1+1) ;
            else
                ind2_start = 1 ;
            end
            
            % --------------------------------------
            % inner loop over images
            for ind2 = ind2_start:N_stacks2
                % ---------------------
                % get filename for second image
                dataFilename_2 = fullfile(bwDir2(ind2).folder, bwDir2(ind2).name) ;
                
                % check if we've already done this calculation
                match_idx_1 = arrayfun(@(x) strcmp(x.label_1,...
                    bwDir1(ind1).name), overlap_struct) ; 
                match_idx_2 = arrayfun(@(x) strcmp(x.label_2,...
                    bwDir2(ind2).name), overlap_struct) ; 
                
                % if we have, move on
                if (sum(match_idx_1 & match_idx_2) > 0)
                    continue
                else
                    % ... otherwise get overlap
                    tic
                    bwMat1 = importdata(dataFilename_1) ;
                    bwMat2 = importdata(dataFilename_2) ;
                    
                    % --------------------
                    % add to structure
                    cc = length(overlap_struct) ;
                    overlap_struct(cc+1).overlap_val = sum(bwMat1(:) & bwMat2(:)) ;
                    overlap_struct(cc+1).label_1 = dataFilename_1 ;
                    overlap_struct(cc+1).label_2 = dataFilename_2 ;
                    toc
                end
            end
            % periodically save so we don't lose data if restart happens
%             if saveFlag
%                 save(overlapFilename, 'overlap_struct')
%             end
        end
        
        % sort results at the end and then save
        labels_comb = arrayfun(@(x) strjoin({x.label_1, x.label_2},'_'),...
            overlap_struct) ; 
        [~, sort_idx] = sort(labels_comb) ; 
        overlap_struct = overlap_struct(sort_idx) ; 
        if saveFlag
            save(overlapFilename, 'overlap_struct')
        end
        overlap_data = overlap_struct ; 
    
    case 'table'
        % -----------------------------------------------------------------
        % check if we already made some progress on overlap calculation
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
                bwDir1 = bwDir1(sort_idx1) ; 
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
                bwDir2 = bwDir2(sort_idx2) ; 
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
            % ---------------------
            % bw image 1
            dataFilename_1 = fullfile(bwDir1(ind1).folder, bwDir1(ind1).name) ;
            lineName1 = label_cell_1{ind1} ; 
            bwMat1 = importdata(dataFilename_1) ;
            
            if strcmp(dataName1,dataName2)
                ind2_start = (ind1+1) ;
            else
                ind2_start = 1 ;
            end
            % make second loop over images
            for ind2 = ind2_start:N_stacks2
                tic
                % ----------------
                % image 2
                dataFilename_2 = fullfile(bwDir2(ind2).folder, bwDir2(ind2).name) ;
                lineName2 = label_cell_2{ind2} ; 
                
                % check to see if this entry is already filled
                if ~isnan(overlap_table{lineName1, lineName2}) && ...
                        ~overWriteFlag
                    fprintf('Already calculated overlap for %s vs %s \n',...
                        lineName1, lineName2)
                    continue
                end
                
                % otherwise load 2nd image
                bwMat2 = importdata(dataFilename_2) ;
                
                % --------------------
                % get overlap
                overlap_table{lineName1, lineName2} = sum(bwMat1(:) & bwMat2(:)) ;
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
        
        % ----------------------------------------------------------------
        %% fill nan entries that shouldn't be there
        % (these would arise if we append a row/column for a new neuron)
        [nan_row, nan_col] = find(isnan(overlap_table{:,:})) ;
        
        % there are lots of intended nan entries for square overlap 
        % matrices (when we compare a set against itself), so make sure to
        % only grab upper triangular elements in those cases
        if strcmp(dataName1,dataName2)
            to_fill_idx = (nan_row < nan_col) ;
            to_fill_row = nan_row(to_fill_idx) ;
            to_fill_col = nan_col(to_fill_idx) ;
        else
            to_fill_row = nan_row ;
            to_fill_col = nan_col ;
        end
        
        % loop through indices that need filling
        for ind = 1:length(to_fill_row)
            tic
            % get appropriate data file for current ROW
            rowNameCurr = rowNamesSort{to_fill_row(ind)} ;
            fn_idx_row = cellfun(@(y) strcmp(y, rowNameCurr), label_cell_1) ;
            dataFilename_row = fullfile(bwDir1(fn_idx_row).folder, ...
                bwDir1(fn_idx_row).name) ;
            
            % get appropriate data file for current COLUMN
            colNameCurr = colNamesSort{to_fill_col(ind)} ;
            fn_idx_col = cellfun(@(y) strcmp(y, colNameCurr), ...
                label_cell_2) ;
            dataFilename_col = fullfile(bwDir2(fn_idx_col).folder, ...
                bwDir2(fn_idx_col).name) ;
            
            % load both data sets
            bwMatRow = importdata(dataFilename_row) ;
            bwMatCol = importdata(dataFilename_col) ;
            
            % get overlap and add to table
            overlap_table{rowNameCurr, colNameCurr} = ...
                sum(bwMatRow(:) & bwMatCol(:)) ;
            toc
            
        end
%         nan_row_ind = find(all_nan_row_idx) ;
%         for indd = nan_row_ind
%             fn_idx = cellfun(@(y) strcmp(
%             dataFilename1 =
%         end
        
        % -------------------------------------------------------------  
        %% save final results?
        if saveFlag
            save(overlapFilename, 'overlap_table')
        end
        
        % give generic name for output
        overlap_data = overlap_table ; 
    otherwise
        fprintf('Invalid data mode selection \n')
        keyboard
end

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