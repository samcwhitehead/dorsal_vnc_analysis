% -------------------------------------------------------------------------
% function to plot matrix + dendrogram for a given set of anatomical
% annotations
% -------------------------------------------------------------------------
function [ax1, ax2, outperm] = plot_connectivity_mat(fig, ax1, ax2, annotationVals, ...
    tree, neuron_labels, target_labels, T_clust, colorThresh, lineColor, ...
    colorClustFlag, xlim, ylim_list, fontSize, fontSizeSmall, ...
    colorPatchFlag, flipYTickLabelFlag, cluster_cmap)
% -----------------------------------
%% inputs and params
if ~exist('T_clust','var') || isempty(T_clust)
    T_clust = ones(size(tree, 1) + 1, 1) ; 
end
if ~exist('colorThresh','var') || isempty(colorThresh)
    colorThresh = Inf ; 
end
if ~exist('lineColor','var') || isempty(lineColor)
    lineColor = [0, 0, 0] ; 
end
if ~exist('colorClustFlag','var') || isempty(colorClustFlag)
    colorClustFlag = false ; 
end
if ~exist('xlim','var') || isempty(xlim)
    xlim = [] ; 
end
if ~exist('ylim_list','var') || isempty(ylim_list)
    ylim_list = [] ; 
end
if ~exist('fontSize','var') || isempty(fontSize)
    fontSize = 6 ; 
end
if ~exist('fontSizeSmall','var') || isempty(fontSizeSmall)
    fontSizeSmall = 4.25 ; 
end
if ~exist('colorPatchFlag','var') || isempty(colorPatchFlag)
    colorPatchFlag = true ; 
end
if ~exist('flipYTickLabelFlag','var') || isempty(flipYTickLabelFlag)
    flipYTickLabelFlag = true ; 
end
if ~exist('cluster_cmap','var') || isempty(cluster_cmap)
    cluster_cmap = 'Dark2' ; 
end

% get axis limits for data
if isempty(xlim)
    xlim = [0.5, length(neuron_labels) + 0.5] ; 
    ylim_list = {[0, max(tree(:,3))], [0.5, length(target_labels) + 0.5]} ; 
end

% other plot params
lw_dend = 0.75 ;

fontName = 'arial' ; 

% ---------------------------------------------------------------
% get colors for plot. 
% for matrix values, this changes depending on whether we're using new or
% old annotations
annot_unique = unique(annotationVals(:)) ; 

% some (darker) colors that we'll use regardless
nothing_color = [1.0, 1.0, 1.0] ; % white
input_color = [33,102,172]/255 ; % blue
output_color = [178,24,43]/255 ; % red 
both_color = [118,42,131]/255 ; % purple

if ~all(round(annot_unique) == annot_unique)
    % new case. need light and dark version of each color for strong/weak
    % inputs/outputs
    weak_input_color = [147,195,255]/255 ; % light blue
    weak_output_color = [255,133,128]/255 ; % light red 
    weak_both_color = [218, 136, 229]/255 ; % light purple
    color_mat = [nothing_color; weak_input_color; input_color; ...
        weak_output_color; output_color; weak_both_color; both_color] ; 
else
    % old case labels are 0 for nothing, 1 for input, 2 for output, and
    % 3 for both ==> white, blue, red, purple
    color_mat = [nothing_color; input_color; output_color; both_color] ;
end
if colorClustFlag
    if ischar(cluster_cmap)
        cluster_colors = brewermap(50,cluster_cmap) ;
    else
        cluster_colors = cluster_cmap ; 
    end
else
    cluster_colors = repmat(lineColor, 50, 1) ; 
end

% patch properties
patch_alpha = 0.3 ; 
% -----------------------------------------
%% plot dendrogram
set(fig,'CurrentAxes',ax1)
hold(ax1,'on')
[h_dend, ~,outperm] = dendrogram( tree,0,'ColorThreshold',colorThresh) ;
set(h_dend, 'Parent',ax1)
close gcf
%set(h_dend,'LineWidth',1.0, 'Color','k') ;
set(ax1,'xlim',xlim, 'ylim', ylim_list{1})
%axis equal

% change line color and width for dendrogram
set(h_dend,'LineWidth', lw_dend) ;
dend_colors = vertcat(h_dend.Color) ;
[dend_colors_unique, iA, iC] = unique(dend_colors, 'rows', 'stable') ;

for q = 1:length(iA)
    color_curr = dend_colors_unique(q,:) ;
    if sum(color_curr) < 0.1
        continue
    end
    change_ind = find(iC == q) ;
    for p = 1:length(change_ind)
        h_dend(change_ind(p)).Color = cluster_colors(q,:) ;
    end
end

% use transparent colored patches to indicate clustering?
if colorPatchFlag
    if ischar(cluster_cmap)
        patch_colors = brewermap(50,cluster_cmap) ; 
    else
        patch_colors = cluster_cmap ; 
    end
    T_clust_unique = unique(T_clust) ;
    N_clust = length(T_clust_unique) ; 
    T_clust_perm = T_clust(outperm) ; 
    
    switch_pts = [0; find(diff(T_clust_perm) ~= 0) ; length(T_clust)] ; 
   patch_array = gobjects(N_clust,1) ; 
   
   y1 = ylim_list{1}(1) ; 
   y2 = ylim_list{1}(2) ; 
   for p = 1:(length(switch_pts) - 1)
       x1 = switch_pts(p) + 0.5 ; 
       x2 = switch_pts(p+1) + 0.5 ; 
       patch_array(p) = patch([x1, x1, x2, x2], [y1, y2, y2, y1], ...
           patch_colors(T_clust_unique(p),:),'EdgeColor','none',...
           'FaceAlpha', patch_alpha) ; 
       
       uistack(patch_array(p), 'bottom')
   end
end
axis(ax1,'off')

%--------------------------------------------------
%% plot matrix image 
ax2 = plot_annotation_mat_toAxis(fig, ax2, annotationVals, ...
    outperm, target_labels, neuron_labels, flipYTickLabelFlag, color_mat,...
    ylim_list{2}, xlim, fontSize, fontSizeSmall, fontName) ; 

end