% -------------------------------------------------------------------------
% function to fetch overlap matrix if we want to compare more than two
% neuron types
% -------------------------------------------------------------------------
function [overlap_table, vox_size_table1, vox_size_table2] = ...
    loadCombinedOverlapMat(dataNames1, dataNames2, rootPath, maskedFlag)
% -------------------------
%% inputs and params
if ~exist('rootPath','var') || isempty(rootPath) 
    [mfilePath, ~, ~] = fileparts(mfilename('fullpath')) ; 
    figDirectory = fileparts(mfilePath) ;
    parentDirectory = fileparts(figDirectory) ; 
    rootPath = fullfile(parentDirectory, 'data') ;
end
if ~exist('maskedFlag','var') || isempty(maskedFlag) 
   maskedFlag = true; 
end
% if ~exist('fillBlocksFlag','var') || isempty(fillBlocksFlag) 
%     % fill in output overlap table so it contains all possible combinations
%     % of elements of dataNames1 and dataNames2?
%    fillBlocksFlag = false; 
% end

% define path to overlap matrices
overlapPath = fullfile(rootPath, 'overlap_calculations') ; 
if maskedFlag
   overlapPath = fullfile(overlapPath,'masked') ;  
end

% general string format for overlap data filenames
fn_str = '%s_%s_overlap_table.mat' ;

% string for vox number table filename and folder
vox_tab_folder = 'binarized_new' ; 
if maskedFlag
   vox_tab_fn = 'vox_size_table_masked.mat' ;
else
   vox_tab_fn = 'vox_size_table.mat' ; 
end
% -----------------------------------------------------------------
%% make sure input data names are in cell format
if ~iscell(dataNames1)
    dataNames1 = {dataNames1} ; 
end
if ~iscell(dataNames2)
    dataNames2 = {dataNames2} ; 
end

% -----------------------------------------------------------------
%% get list of filenames for matrices and load
[xx, yy] = meshgrid(1:length(dataNames1), 1:length(dataNames2)) ;
comb_ind = [xx(:), yy(:)] ;

% number of combinations
N_comb = size(comb_ind,1) ;

% use combinations of neuron types to get possible filenames
filename_list = cell(N_comb,1) ; % initialize storage
overlap_table_list = cell(N_comb,1) ; 

% loop over combinations to get file names and load tables
for k = 1:N_comb
    % current filename
    fn_curr = sprintf(fn_str, dataNames1{comb_ind(k,1)},...
        dataNames2{comb_ind(k,2)}) ;
    
    % add filename to list
    filename_list{k} = fn_curr ; 
    
    % load overlap table
    overlap_table_list{k} = importdata(fullfile(overlapPath, fn_curr)) ; 
end

% -------------------------------------------------------------------
%% try to put overlap tables together
% going to attempt to build by column, which corresponds to dataName2 (this
% is sort of arbitrary, but matches the test case I want to look at)
N_col_blocks = length(unique(comb_ind(:,2))) ;
overlap_table_col_blocks = cell(1,N_col_blocks) ; 

% start by looping over all combinations that involve the same dataName2
for m = 1:N_col_blocks
    % find indices in comb_ind that correspond to current dataName2
    ind_curr = find(comb_ind(:,2) == m) ; 
    
    % load first overlap table 
    overlap_table_curr = overlap_table_list{ind_curr(1)} ; 
    
    % get column names for this first overlap table -- will use these for
    % all
    colNames = overlap_table_curr.Properties.VariableNames ; 
    
    % loop over other combinations with this dataName2 and try to
    % concatenate them on
    for n = 2:length(ind_curr)
        % read out current overlap table
        table_curr = overlap_table_list{ind_curr(n)} ; 
        
        % set column names to standard from first table
        if size(table_curr,2) == size(overlap_table_curr,2)
            table_curr.Properties.VariableNames = colNames ;
%         elseif size(table_curr,2) > size(overlap_table_curr,2)
%             colNamesCurr = table_curr.Properties.VariableNames ; 
%             col_match_idx = ismember(colNamesCurr, colNames) ; 
        else
            keyboard
        end
        % concatenate tables (vertically)
        overlap_table_curr = [overlap_table_curr ; table_curr] ;  
    end
    
    % add resulting column block to cell array 
    overlap_table_col_blocks{m} = overlap_table_curr ; 
end

% --------------------------------------------------------------------
%% concatenate column blocks together 
% NB: this should be one line, but there are a bunch of issues with naming
% inconsistencies
overlap_table = overlap_table_col_blocks{1} ;
rowNames = overlap_table.Properties.RowNames ; % just get a consistent set of row names
for q = 2:N_col_blocks
    % read out current table block
    col_block_curr = overlap_table_col_blocks{q} ; 
    
    % assign it the rowNames we read out before, to keep things consistent
    col_block_curr.Properties.RowNames = rowNames ; 
    
    % add to aggregated table
    overlap_table = [overlap_table, col_block_curr] ;
end

% ----------------------------------------------------------------------
%% need to also get vcx number tables
% for dataName1 list
for ii = 1:length(dataNames1) 
    % get path to current vox num table
    dataPath = fullfile(rootPath, 'vox_size_tables', dataNames1{ii}, ...
        vox_tab_fn) ; 
    
    % load it
    vox_size_table_curr = importdata(dataPath) ; 
    
    % concatenate
    if ii == 1
       vox_size_table1 = vox_size_table_curr ; 
    else
        vox_size_table1 = [vox_size_table1 ; vox_size_table_curr] ; 
    end
end

% for dataName2 list
for jj = 1:length(dataNames2) 
    % get path to current vox num table
    dataPath = fullfile(rootPath, 'vox_size_tables', dataNames2{jj}, ...
        vox_tab_fn) ; 
    
    % load it
    vox_size_table_curr = importdata(dataPath) ; 
    
    % concatenate
    if jj == 1
       vox_size_table2 = vox_size_table_curr ;
    else
        vox_size_table2 = [vox_size_table2 ; vox_size_table_curr] ; 
    end
end

end