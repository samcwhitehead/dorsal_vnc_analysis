% -------------------------------------------------------------------------
% function to reorder matrix rows/cols based on hierarchical clustering
% -------------------------------------------------------------------------
function [mat_out, row_sort_ind, col_sort_ind] = ...
    myMatrixClusterSort(mat_in, cluster_dim, distance_type, linkage_type, ...
    leaf_order_criteria)
% ----------------------------
%% inputs and params
if ~exist('cluster_dim','var') || isempty(cluster_dim)
    cluster_dim = 'all' ; % 'all' | 'row' | 'col'
end
if ~exist('distance_type','var') || isempty(distance_type)
    distance_type = 'correlation' ; % can be any method accepted by pdist
end
if ~exist('linkage_type','var') || isempty(linkage_type)
    linkage_type = 'average' ; % 'average' | 'complete' | ...
end
if ~exist('leaf_order_criteria','var') || isempty(leaf_order_criteria)
    % criteria for optimal leaf order
    leaf_order_criteria = 'adjacent' ;  % 'group'
end

% distance types that cause nan values from pdist for all zero rows
dist_type_check_list = {'correlation','cosine','jaccard','hamming',...
    'spearman'} ; 
% ----------------------------------------------------------------------
%% clean up input matrix
mat_out = mat_in ;

% remove any nan entries
nan_idx = isnan(mat_out) ; 
mat_out(nan_idx) = 0 ; 

% for some distance metrics, rows/cols of all zeros give nan values from
% pdist
if ismember(distance_type,dist_type_check_list)
   skip_row_idx = (sum(abs(mat_out),2) == 0) ; 
   skip_col_idx = (sum(abs(mat_out),1) == 0) ; 
else
    skip_row_idx = false(size(mat_out,1),1) ; 
    skip_col_idx = false(1,size(mat_out,2)) ; 
end

% -----------------------------------------------------------------
%% perfrom sorting (depending on dim input)
switch cluster_dim
    case 'all' 
        % in this case, sort columns and then rows of matrix
        % ---------------------------
        % first sort columns
        % ---------------------------
        % *** note the transpose in the argument of "clusterSortRows"
        col_sort_ind = clusterSortRows(mat_out', distance_type, ...
            linkage_type, leaf_order_criteria, skip_col_idx') ;
        
        % apply sorting to matrix
        mat_out = mat_out(:,col_sort_ind) ; 
        
        % ---------------------------
        % then sort rows
        % ---------------------------
        row_sort_ind = clusterSortRows(mat_out, distance_type, ...
            linkage_type, leaf_order_criteria, skip_row_idx) ;
        
        % apply sorting to matrix
        mat_out = mat_out(row_sort_ind,:) ; 
        
    case 'row'
        % in this case, sort only rows of matrix
        row_sort_ind = clusterSortRows(mat_out, distance_type, ...
            linkage_type, leaf_order_criteria, skip_row_idx) ;
        
        % apply sorting to matrix
        mat_out = mat_out(row_sort_ind,:) ; 
        
        % give column index (no sort)
        col_sort_ind = 1:size(mat_out,2) ; 
        
    case 'col'
        % in this case, sort only columns of matrix
        % *** note the transpose in the argument of "clusterSortRows"
        col_sort_ind = clusterSortRows(mat_out', distance_type, ...
            linkage_type, leaf_order_criteria, skip_col_idx') ;
        
        % apply sorting to matrix
        mat_out = mat_out(:,col_sort_ind) ; 
        
        % give column index (no sort)
        row_sort_ind = 1:size(mat_out,1) ; 
    otherwise
        fprintf('Error: invalid cluster dimension \n') 
        keyboard
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HELPER FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ---------------------------------------------------------------------
%% function to sort rows of a matrix
function row_sort_ind = clusterSortRows(mat, distance_type, linkage_type, ...
    leaf_order_criteria, skip_row_idx)
% only sort non-skip rows
mat_curr = mat(~skip_row_idx,:) ; 

% get distance matrix
D = pdist(mat_curr, distance_type) ; 

% for really sparse overlap matrices, we can't really sort, so just return
% original sorting
if isempty(D)
    row_sort_ind = 1:size(mat,1) ; 
    return
end
% deal with any nan or inf values in distance matrix
bad_idx = isnan(D) | isinf(D) ;
D_max = max(D,[],'all','omitnan') ; 
% NB: giving these "bad" entries a very large distance value, but factor 
% is arbitray. should fix
D(bad_idx) = 100*D_max ; 

% get linkage
tree = linkage(D, linkage_type) ; 


% get leaf order 
leaf_order = optimalleaforder(tree,D,'Criteria',leaf_order_criteria,...
    'Transformation','inverse') ; 

% combine sorted non-skip rows and unsorted skip rows
nonskip_row_ind = find(~skip_row_idx) ;
skip_row_ind = find(skip_row_idx) ;
if ~isempty(skip_row_ind)
    row_sort_ind = [nonskip_row_ind(leaf_order); skip_row_ind] ;
else
    row_sort_ind = leaf_order ; 
end

end