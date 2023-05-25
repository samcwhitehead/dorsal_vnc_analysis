% -------------------------------------------------------------------------
% function to plot anatomical annotation digraph
% -------------------------------------------------------------------------
function [ax, h] = plotAnnotationDigraph(ax, annotationMat, target_labels,...
    edgeColor, layout, binEdgeFlag, edgeTrimThresh, legendPts, rmNodesFlag, ...
    normEdgeFlag, removeNeuropils)
% ------------------------------------------------------------
%% inputs and params
% color for edges
if ~exist('edgeColor', 'var') || isempty(edgeColor)
    edgeColor = [0, 0.447, 0.741] ; 
end
% layout of graph
if ~exist('layout', 'var') || isempty(layout)
    layout = 'custom' ; 
end
% make every neuron an edge, or combine them?
if ~exist('binEdgeFlag','var') || isempty(binEdgeFlag)
   binEdgeFlag = true ; 
end
% smallest allowed width for edges--smaller ones will be trimmed. if empty,
% no trimming occurs
if ~exist('edgeTrimThresh','var') || isempty(edgeTrimThresh)
   edgeTrimThresh = [] ; 
end
% two points (neuron numbers) that we'll use to make a legend for edge
% thickness. if empty, don't make legend
if ~exist('legendPts','var') || isempty(legendPts)
   legendPts = [] ;  
end
% remove nodes with no edges?
if ~exist('rmNodesFlag','var') || isempty( rmNodesFlag)
    rmNodesFlag = false ;  
end
% normalize edge width to max current degree?
if ~exist('normEdgeFlag','var') || isempty(normEdgeFlag)
    normEdgeFlag = true ;  
end
% remove ANm, VAC, and mVAC nodes to simplify graph?
if ~exist('removeNeuropils','var') || isempty(removeNeuropils)
    removeNeuropils = {'ANm','VAC','mVAC','Y'} ;  % cell(0) ; 
end 

% --------------------------
% graph node properties

% if strcmp(layout,'custom')
%     nodeSize = 3 ; 
% else
%     nodeSize = 5 ;  % 5
% end
nodeSize = 3 ;  % 5
edgeScale = 5 ; % edge weight is divided by this! % 4
nodeColor = 0.3*[1 1 1] ; 
nodeFontSize = 6 ; 
edgeAlpha = 0.5 ; 

outlineColor = 0.*[1,1,1] ; 
outlineWidth = 0.75 ; 

arrowSize = 6 ; 
arrowPosition = 0.75 ; 
% ---------------------------------------------------
% set of rare neuropils that we might need to remove depending on
% rmRareNeuropilsFlag


% -------------------------------------------------------------
%% get positions of target vnc regions
% remove nodes like ANm, VAC, and mVAC?
if ~isempty(removeNeuropils)
    % find where neuropils we want to remove show up
    [~, rm_np_ind] = ismember(removeNeuropils, target_labels) ; 
    rm_np_ind = rm_np_ind(rm_np_ind~=0) ; % remove indices where we don't have matches
    rm_np_idx = false(size(target_labels)) ;  % convert to logical idx
    rm_np_idx(rm_np_ind) = true ; 
    
    % remove these elements from target_label list and annotation mat
    target_labels = target_labels(~rm_np_idx) ; 
    annotationMat = annotationMat(:, ~rm_np_idx) ; 
end

% for VNC layout
vnc_template_struct = get_vnc_layout_from_im(target_labels) ;
target_positions = vnc_template_struct.neuropil_locs ; 

% for circle layout (fixed positions of neuropil regions in circular
% arrangement)
circle_x = [1.0000, 0.8660, 0.5000, 0, -0.5000, -0.8660, -1.0000, ...
    -0.8660, -0.5000, 0, 0.5000, 0.8660] ;
circle_y = [0, 0.5000, 0.8660, 1.0000, 0.8660, 0.5000, 0, -0.5000,...
    -0.8660, -1.0000, -0.8660, -0.5000] ;

% ------------------------------------------------------------
%% convert binary matrix into source and target node pairs
s = [] ; 
t = [] ; 
for i = 1:size(annotationMat,1)
    neuronAnnCurr = annotationMat(i,:) ; 
    inputs = find((neuronAnnCurr == 1) | (neuronAnnCurr == 3)) ; 
    outputs = find((neuronAnnCurr == 2) | (neuronAnnCurr == 3)) ; 
    
    % find pairs
    [idx2, idx1] = find(true(numel(outputs),numel(inputs)));
    s = [s ; reshape(inputs(idx1), [], 1)] ; 
    t = [t ; reshape(outputs(idx2), [], 1)] ;   
end

% ------------------------------------------------------------
%% if we're binning, do it now
if binEdgeFlag
   st_mat = [s, t] ; 
   [st_unique, ~, ic] = unique(st_mat,'rows','stable') ; 
   % Count Occurrences
   weights = accumarray(ic, 1);                            
   s = st_unique(:,1) ;
   t = st_unique(:,2) ; 
else
    weights = ones(size(s)) ; 
end

% -------------------------------------------------------------
%% make graph and plot
% get edge widths based on weights 
G = digraph(s,t,weights,target_labels) ; 

edge_weights = G.Edges.Weight ; 

% normalize to max?
if normEdgeFlag
    edge_weights_norm = ...
        (edge_weights - min(edge_weights))./max(edge_weights) ;
    edge_widths = 4*edge_weights_norm + 0.1 ;
    
    % -----------------------------
    % trim small edges?
    if ~isempty(edgeTrimThresh)
        small_edge_idx = (edge_widths < edgeTrimThresh) ;
        G = rmedge(G, find(small_edge_idx)) ;
        edge_widths = edge_widths(~small_edge_idx) ;
        edge_weights = G.Edges.Weight ;
    end
else
    % if we're not normalizing, just take SCALED weight as width
    edge_widths = edge_weights./edgeScale ;
end

% ---------------------------------------------
% remove nodes with no edges?
if rmNodesFlag
    % get number of edges in and out of each node
   in_deg = indegree(G) ; 
   out_deg = outdegree(G) ; 
   
   % find nodes with zero in or out edges
   rm_idx = (in_deg < 1) & (out_deg < 1) ;
   
   % remove these nodes
   G = rmnode(G, find(rm_idx)) ; 
   
   % if using a graph layout where we set node position based on fixed
   % coordinates, remove corresponding coordinates
   if length(circle_x) == length(rm_idx)
       circle_x = circle_x(~rm_idx) ;
       circle_y = circle_y(~rm_idx) ;
   end
   if size(target_positions,1) == length(rm_idx)
       target_positions = target_positions(~rm_idx,:) ;
   end
end

% --------------------------------------------------
% make plot based on layout selection
switch layout
    case 'custom'
        hold(ax,'on')
        vnc_outline = vnc_template_struct.vnc_outline ; 
        plot(vnc_outline(:,1), vnc_outline(:,2), '-',...
            'LineWidth', outlineWidth, 'Color', outlineColor) 
        
        h = plot(ax, G, 'LineWidth',edge_widths, 'nodeColor', nodeColor, ...
            'markerSize',nodeSize, 'EdgeColor', edgeColor, ...
            'Xdata', target_positions(:,1), 'Ydata', target_positions(:,2),...
            'NodeFontSize',nodeFontSize, 'NodeFontName','arial') ;
       
        
    case 'custom circle'
        corner_x = 1.5 ; % 2.0 
        corner_y = -1.0 ; 
        x_data = [circle_x, corner_x] ; 
        y_data = [circle_y, corner_y] ; 
        h = plot(ax, G, 'LineWidth',edge_widths, 'nodeColor', nodeColor, ...
            'markerSize',nodeSize, 'EdgeColor', edgeColor, ...
            'Xdata', x_data, 'Ydata', y_data,...
            'NodeFontSize',nodeFontSize, 'NodeFontName','arial') ;
    case 'circle and patch'
        h = plot(ax, G, 'LineWidth',edge_widths, 'nodeColor', nodeColor, ...
            'markerSize',nodeSize, 'EdgeColor', edgeColor, ...
            'Layout','circle', 'NodeFontSize',nodeFontSize, ...
            'NodeFontName','arial') ;
        h.NodeLabelMode = 'manual' ; 
%         h.XData = circle_x ; 
%         h.YData = circle_y ; 
%         h = plot(ax, G, 'LineWidth',edge_widths, 'nodeColor', nodeColor, ...
%             'markerSize',nodeSize, 'EdgeColor', edgeColor, ...
%             'Xdata', circle_x, 'Ydata', circle_y,...
%             'NodeFontSize',nodeFontSize, 'NodeFontName','arial') ;
    case 'circle'
        h = plot(ax, G, 'LineWidth',edge_widths, 'nodeColor', nodeColor, ...
            'markerSize',nodeSize, 'EdgeColor', edgeColor, ...
            'Layout',layout, 'NodeFontSize',nodeFontSize, 'NodeFontName','arial') ;
    case 'force'
        h = plot(ax, G, 'LineWidth',edge_widths, 'nodeColor', nodeColor, ...
            'markerSize',nodeSize, 'EdgeColor', edgeColor, ...
            'Layout',layout, 'NodeFontSize',nodeFontSize, 'NodeFontName','arial')  ;
    otherwise
        fprintf('Invalid layout selection: %s \n',layout)
        keyboard
end

% set arrow size/position
h.ArrowSize = arrowSize ; 
h.ArrowPosition = arrowPosition ; 

% set edge alpha
h.EdgeAlpha = edgeAlpha ; 

% TEMP -- try to use linestyle to make some edges less visible
for k = 1:length(edge_weights)
    if edge_weights(k) <= 2
        highlight(h, G.Edges.EndNodes(k,:), 'LineStyle',':') ;
    end
end

% TEMP testing removing marker fill color for nods w/o input or output
in_deg = indegree(G) ;
out_deg = outdegree(G) ;
no_io_idx = (in_deg < 1) & (out_deg < 1) ;
no_io_ind = find(no_io_idx) ;
highlight(h, no_io_ind, 'Marker','none') ;


% make non-VNC node a little different
non_vnc_ind = find(cellfun(@(y) strcmp(y, 'nonVNC'), target_labels)) ; 
if numel(non_vnc_ind) > 0
    highlight(h,target_labels{non_vnc_ind}, 'NodeColor','k')
end

axis(ax, 'equal')
axis(ax ,'off')
set(ax,'FontSize', nodeFontSize, 'FontName', 'arial') 
% ----------------------------------
% make legend?
if ~isempty(legendPts)
    % make cell array for string versions of legendPts
    legendPtLabels = arrayfun(@(x) [num2str(x) '  neurons'], legendPts,...
        'UniformOutput',0) ; 
    % make interpolant to fit legend points to data points
    edge_weights_unique = unique(edge_weights) ; 
    edge_widths_unique = unique(edge_widths) ; 
    line_widths = interp1(edge_weights_unique, edge_widths_unique, ...
        legendPts,'linear','extrap') ; 
    
    % plot dummy variables to make custom legend
    N_legend_pts = length(line_widths) ;
    h_array = gobjects(N_legend_pts,1) ; 
    
    hold(ax, 'on')
    for k = 1:N_legend_pts
        h_array(k) = plot(ax, [nan, nan], [nan, nan], '-', ...
            'Color', [edgeColor, edgeAlpha],...
            'LineWidth', line_widths(k)) ;
    end
    
    % make legend
    legend(h_array, legendPtLabels, 'location', 'southwest')   
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HELPER FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------------------------------------------------------
%% function to get VNC outline and neuropil centroids from image
% to be used for making drigraphs that imitate vnc layout
function vnc_template_struct = get_vnc_layout_from_im( neuropil_labels, ...
    imPath, saveFlag, loadFlag)
% ----------------------------------------------------
% inputs and params
% ----------------------------------------------------
if ~exist('neuropil_labels', 'var') || isempty(neuropil_labels)
    % names of neuropil locations in our graph
    neuropil_labels = {'NTct', 'WTct', 'HTct', 'IntTct', 'LTct', 'T1Leg', ...
        'T2Leg', 'T3Leg','AMNp', 'ANm', 'VAC', 'mVAC'} ;
end
if ~exist('imPath', 'var') || isempty(imPath)
    % path where we can find vnc template image. should be in same
    % directory as function
    path_full = mfilename('fullpath') ; 
    [imPath, ~, ~] = fileparts(path_full) ; 
end
if ~exist('saveFlag', 'var') || isempty(saveFlag)
    % save the results of this calculation so we don't have to re-do every
    % time we make a digraph?
    saveFlag = true ; 
end
if ~exist('loadFlag', 'var') || isempty(loadFlag)
    % given the chance, should we load previously calculated results?
    loadFlag = false ; 
end

% number of neuropils we're looking for
N_np = length(neuropil_labels) ; 

% filename for template image and to save output to
% fn = 'digraph_vnc_layout_template_adjusted.png' ; 
fn = 'digraph_vnc_layout_template_new_annotations.png' ;
imPathFull = fullfile(imPath, fn) ;
save_fn = 'vnc_template_struct_new_annot.mat' ; 

% -------------------------------------------------------------
% check if we've already done calculation -- if so, load it?
% -------------------------------------------------------------
if exist(fullfile(imPath,save_fn),'file') && loadFlag
   vnc_template_struct = importdata( fullfile(imPath,save_fn)) ; 
   return
end

% ----------------------------------------------------
% load template image
% ----------------------------------------------------
vnc_template = imread(imPathFull) ; 

% convert to grayscale/8bit
vnc_template = uint8(rgb2gray(vnc_template)) ; 

% ----------------------------------------------------
% perform multithresholding to get neuropil regions
% ----------------------------------------------------
% NB: multithresh doesn't give perfect results because colors aren't evenly
% spaced, so just going to use histcounts
edges = -0.5:255.5 ; 
centers = (edges(1:end-1) + edges(2:end))./2 ;
vnc_template_flat = vnc_template(:) ; 
[N_pix, ~, bins] = histcounts(vnc_template_flat, edges) ; 

% due to something with the loading/processing, some pixel values get
% smeared, so take only bins that correspond to the neuropils and outline
[~, sort_ind] = sort(N_pix, 'descend') ; 
centers_sort = centers(sort_ind) ; 

% get bins for neurpils + outline (skip over background)
bins_keep = centers_sort(2:(N_np + 2)) + 1; % need to add one due to bin index offset

% re-sort bins that we're keeping so we stay in neuropil order
bins_keep = sort(bins_keep, 'ascend') ; 

% -----------------------------------------------------------------
% generate segmented image
% -----------------------------------------------------------------
% initialize label image
vnc_template_label = zeros(size(vnc_template)) ; 

% loop over bins (i.e. regions we want to keep except outline) and add to 
% label image
for k = 1:length(bins_keep)
    % index for current object
   idx_curr = (bins == bins_keep(k)) ; 
   
   % convert index into binary image so we can clean up
   bw_curr = false(size(vnc_template)) ; 
   bw_curr(idx_curr) = true ; 
   bw_curr = bwareaopen(bw_curr, 200) ; 
   
   % add to label image
   vnc_template_label = vnc_template_label + k.*double(bw_curr) ; 
end

% -----------------------------------------------------------------
% get each object in label image
% -----------------------------------------------------------------
% use regionprops to get info about each region
label_props = regionprops(vnc_template_label, 'Centroid') ; 

% for everything but outline, get centroid coordinates (we'll use these for
% node locations)
neuropil_locs = vertcat(label_props(2:end).Centroid);

% for outline, get pixel coordinates
vnc_outline_bw = false(size(vnc_template_label)) ; 
vnc_outline_bw(vnc_template_label == 1) = true ; 
vnc_bw_fill = imfill(vnc_outline_bw,'holes') ;
vnc_outline_perim = bwperim(vnc_bw_fill) ;
%vnc_outline_perim = bwmorph(vnc_outline_perim,'thin', Inf) ; 
vnc_outline_idx = find(vnc_outline_perim == 1) ; 
[vnc_y, vnc_x] = ind2sub(size(vnc_template_label), vnc_outline_idx) ;
vnc_outline = [vnc_x, vnc_y] ; 

% need to sort outline coordinates so we can plot as a line. do this by
% iterative distance search
D = pdist(vnc_outline) ; % get distance between all points
Z = squareform(D) ;

% set diagonal of distance matrix to infinity so it doesn't throw off search
N = size(vnc_outline,1) ; % number of pixels
Z(1:N+1:end) = Inf ; 

% initialize storage for sort index
idx = 1 ; 
out_idx = zeros(N,1) ; 
out_idx(1) = idx ; 

% loop through and find closest pixel
for m = 2:N
    % start ind is the point we're starting at -- want to find closest
    % point to it
    start_ind = idx ; 
    
    % find index of minimum distance point 
    [~,idx] = min(Z(start_ind,:));
    
    % kill off the column corresponding to previous closest point so we
    % don't pass it again
    Z(:,start_ind) = Inf;
    % update search row to the one we just found
    out_idx(m) = idx;
end

% sort vnc outline with new order
vnc_outline = vnc_outline(out_idx,:) ; 

% check if we still have any large jumps
vnc_outline_diff = [0; myNorm(diff(vnc_outline))]; 
jump_idx = (vnc_outline_diff > sqrt(2)) ; 
vnc_outline(jump_idx,:) = nan ; 

if (0)
   figure ; 
   RGB = label2rgb(vnc_template_label) ; 
   imshow(RGB) 
   hold on
   plot(neuropil_locs(:,1), neuropil_locs(:,2), 'ko', 'MarkerFaceColor','k')
   plot(vnc_outline(:,1), vnc_outline(:,2), 'k-','LineWidth',2)
end

% ------------------------------------------------------
% convert to cartesian coordinates and add to struct
% ------------------------------------------------------
% coordinate conversion
neuropil_locs(:,2) = size(vnc_template,1) - neuropil_locs(:,2) + 1;
vnc_outline(:,2) = size(vnc_template,1) - vnc_outline(:,2) + 1;

% add to struct
vnc_template_struct = struct() ; 
vnc_template_struct.neuropil_locs = neuropil_locs ; 
vnc_template_struct.vnc_outline = vnc_outline ; 


% -------------------------------------
% save results?
% -------------------------------------
if saveFlag
   savePathFull = fullfile(imPath, save_fn) ; 
   save(savePathFull, 'vnc_template_struct') ; 
end
end

%{
target_positions = [ 227, 112 ; ... % neck neuropil
                     161, 110 ; ... % wing neuropil
                     95 , 111 ; ... % haltere "
                     153, 97 ; ... % tectulum
                     174, 87 ; ... % lower tectulum
                     239, 81 ; ... % T1 leg
                     155, 71 ; ... % T2 leg
                     58 , 80 ; ... % T3 leg
                     197, 75 ; ... % AMN
                     17 , 100 ; ... % AS
                     154, 52 ; ... % VAC
                     186, 66 ] ; % ... % mVAC
%                      161, 125 ] ;  % nonVNC
%}