% -------------------------------------------------------------------------
% function to cluster neurons by their input/output patterns. this was
% previously just done in the plotting function, but I think should be in
% separated
% -------------------------------------------------------------------------
function [conn_clust_struct, conn_clust_table] = ...
    cluster_io_annotations(annotationVals, neuron_labels, ID_nums, ...
    target_labels, max_k, distanceType, linkageMethod, debugFlag, ...
    capitalFactor, useIpsiContraFlag, rmWeakAnnotFlag)
% --------------------
%% inputs and params
if ~exist('max_k','var') || isempty(max_k)
    max_k = 16 ; % maximum number of clusters to evaluate
end
if ~exist('distanceType','var') || isempty(distanceType)
    distanceType = 'correlation'  ; % distance metric to use when comparing i/o patterns
end
if ~exist('linkageMethod','var') || isempty(linkageMethod)
    linkageMethod = 'complete' ; % linkage rule for dendrogram (hierarchical clust)
end
if ~exist('debugFlag','var') || isempty(debugFlag)
    debugFlag = false ; % plot some diagnostics?
end
if ~exist('capitalFactor','var') || isempty(capitalFactor)
    % FOR NEW ANNOTATIONS -- how much to weight upper vs lowercase letters
    capitalFactor = 2 ; 
end
if ~exist('useIpsiContraFlag','var') || isempty(useIpsiContraFlag)
    % FOR NEW ANNOTATIONS -- keep ipsi and contra annotations separate?
    useIpsiContraFlag = true; 
end
if ~exist('rmWeakAnnotFlag','var') || isempty(rmWeakAnnotFlag)
    % FOR NEW ANNOTATIONS -- only keep strong (upper case) annotations?
    rmWeakAnnotFlag = true; 
end

evaluationMethod = 'Gap' ; % 'Gap' % used to determine optimal cluster number
% ------------------------------------------------------
%% reshape annotation values for distance calculations
if iscell(annotationVals)
    % convert character pairs to numbers. annotationValsMat will store the
    % NxMx(2 or 4) array, not reshaped
    [annotationValsReshape, annotationValsMat] = ...
        convertNewAnnotationVals(annotationVals, capitalFactor, ...
        useIpsiContraFlag, rmWeakAnnotFlag) ; 
    
else
    % turn matrix of 1,2,3 into a larger, binary matrix
    annotationValsReshape = reshapeAnnotationMat(annotationVals) ;
    % need a new variable for storing annotation vals to add to struct
    annotationValsMat = annotationVals ; 
end

% ---------------------------------------------------------------------
%% try to determine optimal cluster number
eva = evalclusters(annotationValsReshape,'linkage', evaluationMethod, ...
    'KList', [3:max_k],'Distance','correlation') ; %,...
    %'ReferenceDistribution','uniform') ;
max_clusts = eva.OptimalK ;

if debugFlag
   figure ; 
   hold on
   plot(eva)
   
   axis tight
   ylim_curr = get(gca,'ylim') ; 
   plot(max_clusts.*[1, 1], ylim_curr, 'r--')
   
   xlabel('Number of clusters')
   ylabel('Criterion Value') 
end
% ---------------------------------------------------------------------
%% perform hierarchical clustering
D = pdist(annotationValsReshape, distanceType) ;
tree = linkage(D,linkageMethod) ;

T_clust = cluster(tree, 'maxclust',max_clusts) ;

% -----------------------------------------------------
%% organize cluster results into struct and table
conn_clust_struct = struct() ;
for i = 1:length(T_clust)
    conn_clust_struct(i).neuron_label = neuron_labels{i} ;
    conn_clust_struct(i).ID = ID_nums{i} ;
    conn_clust_struct(i).cluster_idx = T_clust(i) ;
    for j = 1:length(target_labels)
        fieldName = ['target_' target_labels{j}] ; 
        conn_clust_struct(i).(fieldName) = squeeze(annotationValsMat(i,j,:)) ;
    end
end

% generate table with the same results
conn_clust_table = struct2table(conn_clust_struct) ;

% add some other misc. stuff we may need (this wouldn't fit in
% table)
conn_clust_struct(1).tree = tree ;
conn_clust_struct(1).linkageMethod = linkageMethod ;        
conn_clust_struct(1).distanceType = distanceType ; 
conn_clust_struct(1).max_k = max_k ; 
end
