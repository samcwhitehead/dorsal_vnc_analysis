% -------------------------------------------------------------------------
% script to combine all overlap tables into one large matrix with labels,
% which can then be exported to a file format readable by cytoscape
% (software for making directed graphs)
%
% NB: just going to assume we're using masked overlap -- if this switches,
% will need to alter code
%
% *** need to swap row/col dimensions (I think) to match cytoscape format
% -------------------------------------------------------------------------
%% path info and params
[mfilePath, ~, ~] = fileparts(mfilename('fullpath')) ;
figDirectory = fileparts(mfilePath) ;
parentDirectory = fileparts(figDirectory) ;
rootPath = fullfile(parentDirectory, 'data') ;
dataPath = fullfile(rootPath, 'processed_anatomy_data') ; 
overlapPath = fullfile(rootPath, 'overlap_calculations', 'masked') ;

% all types of neurons used in overlap analyses
neuronTypes = {'IN', 'DN', 'VUM', 'MN'} ;

% save output?
saveFlag = true ;

% where to save output
% savePath = fullfile(overlapPath, 'full_graph') ;
savePath = fullfile(overlapPath, 'combined') ;
if ~exist(savePath, 'dir') && saveFlag
    mkdir(savePath)
end

transposeFlag = false ;

% --------------------------------------------------------
%% enumerate possible overlap combinations
[xx, yy] = meshgrid(1:length(neuronTypes), 1:length(neuronTypes)) ;
comb_ind = [xx(:), yy(:)] ;

% number of combinations
N_comb = size(comb_ind,1) ;

% use combinations of neuron types to get possible filenames
filename_list = cell(N_comb,1) ; % initialize storage
fn_str = '%s_%s_overlap_table.mat' ;

% loop over combinations
for k = 1:N_comb
    filename_list{k} = sprintf(fn_str, neuronTypes{comb_ind(k,1)},...
        neuronTypes{comb_ind(k,2)}) ;
end

% ------------------------------------------------------------------
%% loop over overlap files and load
% where we don't have an overlap file, we'll assume it's all zeros
% going to try to do this by concatenating tables first -- hopefully that
% will make keeping track of the labels less of a headache
% ------------------------------------------------------------------
% initialize counter
cc = 1 ;

% double loop over neuron types
for ii = 1:length(neuronTypes)
    % current row neurons
    neuronType1 = neuronTypes{ii} ;
    
    for jj = 1:length(neuronTypes)
        tic
        % current column neurons
        neuronType2 = neuronTypes{jj} ;
        % get full file path for current overlap
        fn_curr = filename_list{cc} ;
        filePathFull = fullfile(overlapPath, fn_curr) ;
        
        % if overlap file exists, load it.
        if exist(filePathFull,'file')
            overlap_table_curr = importdata(filePathFull) ;
        else
            % otherwise, we need to generate a table of zeros with the
            % appropriate size and labels. To do that, first get data
            % directories for neuron type:
            dataDir1 = dir(fullfile(dataPath,neuronType1,'binarized_sym',...
                '*_coords.mat'));
            dataDir2 = dir(fullfile(dataPath,neuronType2,'binarized_sym',...
                '*_coords.mat'));
            
            % how many data files in each directory?
            N_stacks1 = length(dataDir1) ;
            N_stacks2 = length(dataDir2) ;
            
            % get labels from each data set
            label_cell_1 = cell(N_stacks1,1) ;
            for k = 1:N_stacks1
                fn = dataDir1(k).name ;
                label_cell_1{k} = getLineName(fn, neuronType1, true) ;
            end
            
            label_cell_2 = cell(N_stacks2,1) ;
            for k = 1:N_stacks2
                fn = dataDir2(k).name ;
                label_cell_2{k} = getLineName(fn, neuronType2, true) ;
            end
            
            % make empty table
            overlap_mat_curr = zeros(N_stacks1,N_stacks2) ; 
            overlap_table_curr = array2table(overlap_mat_curr) ; 
            
            % assign labels as row and variable names
            overlap_table_curr.Properties.RowNames = label_cell_1 ; 
            overlap_table_curr.Properties.VariableNames = label_cell_2 ; 
        end
        
        % ---------------------------------------------------------------
        % we're constructing full matrix in block rows, so if this is the
        % first entry for this row, make it the current table. otherwise
        % concatenate. First need to sort rows and columns so that
        % everything matches up consistently though
        rowNames = overlap_table_curr.Properties.RowNames ; 
        colNames = overlap_table_curr.Properties.VariableNames ; 
        
        % * seem to be consistently dealing with issue of underscore vs
        % hyphen in row/col names. going to try a hacky fix for that:
        if strcmp(neuronType1,'IN')
            % loop over ROWS
            for kk = 1:length(rowNames)
               rowNameCurr = rowNames{kk} ; 
               rowNameSplit = strsplit(rowNameCurr,'_') ; 
               if length(rowNameSplit) > 2 
                  rowNameNew = strjoin([strjoin(rowNameSplit(1:2),'-'),...
                      rowNameSplit(3:end)],'_') ;
                  rowNames{kk} = rowNameNew ; 
               end
            end
        end
        if strcmp(neuronType2,'IN')
            % loop over COLUMNS
            for kk = 1:length(colNames)
               colNameCurr = colNames{kk} ; 
               colNameSplit = strsplit(colNameCurr,'_') ; 
               if length(colNameSplit) > 2 
                  colNameNew = strjoin([strjoin(colNameSplit(1:2),'-'),...
                      colNameSplit(3:end)],'_') ;
                  colNames{kk} = colNameNew ; 
               end
            end
        end
        
        % put new row names back into table
        overlap_table_curr.Properties.RowNames = rowNames ; 
        overlap_table_curr.Properties.VariableNames = colNames ; 
        
        % sort row and column names
        row_sort_ind = sort(rowNames) ;
        col_sort_ind = sort(colNames) ; 
        
        overlap_table_curr = overlap_table_curr(row_sort_ind, col_sort_ind) ; 
        
        % add neuron type as variable description (columns)
        colNeuronTypes = cell(size(overlap_table_curr,2),1) ; 
        colNeuronTypes(:) = {neuronType2} ; 
        overlap_table_curr.Properties.VariableDescriptions = colNeuronTypes ;
        
        % intialize or add to current table block row
        if jj == 1
            overlap_table_row = overlap_table_curr ;
        else
            try
                overlap_table_row = [overlap_table_row, overlap_table_curr] ;
            catch
                fprintf('Error concatenating! \n')
                keyboard
            end
        end
        
        % print update
        toc
        fprintf('Completed %s (%d/%d) \n', fn_curr, cc, N_comb)
        
        % increment counter 
        cc = cc + 1 ; 
    end
    
    
    % after completing a block row: if it's the first, use it to initialize
    % full output. otherwise concatenate. 
    if ii == 1
        % create full overlap table
        overlap_table_all = overlap_table_row ;

    else
        % tag on current block row
        try
            overlap_table_all = [overlap_table_all ; overlap_table_row] ;
        catch
            keyboard
        end
        
    end
      
end

% --------------------------------------------
%% take table transpose?
if transposeFlag
   % not sure what will happen to variable descriptions if I use 
   % myTableTranspose, so just going to take transpose of values (row and
   % col labels should match)
   overlap_table_all{:,:} = overlap_table_all{:,:}' ; 
end
% -------------------------------------------
%% save output?
if saveFlag
   savePathFull = fullfile(savePath, 'overlap_table_all.mat') ; 
   save(savePathFull,'overlap_table_all')
end
