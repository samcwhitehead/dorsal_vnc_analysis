% -------------------------------------------------------------------------
% function to normalize matrix of anatomical overlap when using neuropil
% regions as mask.
%
% could have probably made the old "normalize_overlap_data." compatible
% with this, but sick of twisting in loops for old code
% -------------------------------------------------------------------------
function overlap_data_norm = normalize_masked_overlap_data(overlap_data,...
    vox_num_data1, vox_num_data2, normType)
% ---------------------
%% inputs
if ~exist('vox_num_data1','var') || isempty(vox_num_data1)
    vox_num_data1 = [] ; % normally a vector conainting voxel count for each stack
    normType = 'max' ;
end
if ~exist('vox_num_data2','var') || isempty(vox_num_data2)
    vox_num_data2 = [] ; % normally a vector conainting voxel count for each stack
    normType = 'max' ;
end
if ~exist('normType','var') || isempty(normType)
    normType = 'max' ; % method of normalization. can be 'max', 'sym_size', 'asym_size'
end

INPUT = 1 ; 
OUTPUT = 2 ; 
% --------------------------------------------------------
%% remove any deleted entries from overlap table/matrix
% get row and column names
rowNames = overlap_data.Properties.RowNames ;
colNames = overlap_data.Properties.VariableNames ;

% find any rows or columns with "deleted" in the name
deleted_row_idx = cellfun(@(y) contains(y, 'delete', 'IgnoreCase',1),...
    rowNames) ;
deleted_col_idx = cellfun(@(y) contains(y, 'delete', 'IgnoreCase',1),...
    colNames) ;

% keep only rows/cols without the "deleted" label
overlap_data = overlap_data(~deleted_row_idx, ~deleted_col_idx) ;


% ---------------------------------
%% initialize normalized output
[N_stacks1, N_stacks2] = size(overlap_data) ;
% switch numbers on voxel size vector? I think my overlap matrices might
% not have consistent order
if (size(vox_num_data1,1) == N_stacks2) && ...
        (size(vox_num_data2,1) == N_stacks1)
    temp = vox_num_data1 ;
    vox_num_data1 = vox_num_data2 ;
    vox_num_data2 = temp ;
end

% initialize output
overlap_norm_init = nan(size(overlap_data)) ;
overlap_data_norm = array2table(overlap_norm_init) ;
overlap_data_norm.Properties.RowNames = ...
    overlap_data.Properties.RowNames ;
overlap_data_norm.Properties.VariableNames = ...
    overlap_data.Properties.VariableNames ;

% --------------------------------------------------------
%% deal with overlap mat diagonal 
% it's more complicated with masked data, and not clear that it couldn't be
% interesting with e.g. a bilateral pair of neurons, but not going to worry
% about it for now
if (N_stacks1 == N_stacks2)
    for i = 1:N_stacks1
        overlap_data{i,i} = nan ;
    end
end
    
% ------------------------------------------------------
%% allow for a few different types of normalization
switch normType
    case 'max'
        % normalize so maximum overlap is 1. does not take into account
        % overall size of expression pattern
        max_overlap = nanmax(overlap_data{:,:},[],'all') ;
        overlap_data_norm{:,:} = overlap_data{:,:}./max_overlap ;
        
    case {'sym_size', 'asym_size'}
        % NB: masked data is inherently asymmetric, so these will both
        % ultimately be the same thing. however, for backwards
        % compatibility, leaving both as an option
        
        % divide the voxel overlap number by the average of the voxel
        % counts for the two stacks being compared. this will lead to a
        % symmetric overlap matrix
        % -----------------------------------------
        % loop over rows
        for r = 1:N_stacks1
            % -------------------
            % vox num for rows:
            % -------------------
            % pull out matching entry
            vox_size1 = vox_num_data1.voxCount(r,INPUT) ;
           
            % ----------------------------------------------
            % loop over columns
            for c = 1:N_stacks2
                % ------------------
                % vox num for cols:
                % ------------------
                vox_size2 = vox_num_data2.voxCount(c, OUTPUT) ;
               
                % ------------------------------------------------------
                % use these results to get mean voxel size and use to
                % normalize
                mean_size = nanmean([vox_size1, vox_size2]) ;
%                 mean_size = max([vox_size1, vox_size2]) ;
                overlap_data_norm{r, c} = overlap_data{r, c}./mean_size ;
            end
        end
        
    case 'none'
        % do nothing
        overlap_data_norm{:,:} = overlap_data{:,:} ;
    otherwise
        disp('Invalid normalization type')
        overlap_data_norm = [] ;
        return
end

end