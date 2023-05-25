% -------------------------------------------------------------------------
% function to get permutation of neurons by how they're placed in the
% dendrogram -- will help with cluster number identification
% -------------------------------------------------------------------------
function [T_clust_perm, T_clust_perm_unique, outperm] = ...
    get_cluster_outperm(conn_clust_struct, colorThresh)
% inputs and params (for dendrogram plot)
if ~exist('colorThresh','var') || isempty(colorThresh)
    colorThresh = Inf ; % no color cutting
end
% number of leaf nodes allowed in plot -- set to 0 to allow max
N_leaf_nodes = 0 ; 

% read tree structure from connectivity struct 
tree = conn_clust_struct(1).tree ; 

% make temporary dendrogram plot to get the right permuation of elements
h_temp = figure ; 
[~, ~,outperm] = dendrogram(tree, N_leaf_nodes,...
    'ColorThreshold',colorThresh) ;

% close temporary plot
close(h_temp)

% use dendrogram permutation to re-order cluster labels 
T_clust = [conn_clust_struct.cluster_idx] ;
T_clust_perm = T_clust(outperm) ; 

% also get unique elements of permuted cluster indices (in order)
T_clust_perm_unique = unique(T_clust_perm, 'stable') ; 
end