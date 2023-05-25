% -------------------------------------------------------------------------
% function to summarize hemilineage identity of interneurons by cluster in
% plot form
% -------------------------------------------------------------------------
function ax = plotClusterHemilineageToAxis(ax, conn_clust_struct, ...
    clusterNum, useSubTypeFlag, titleFlag, colorMode, scaleChartFlag, ...
    sortWedgeFlag)
% ---------------------
%% inputs and params
if ~exist('ax','var') || isempty(ax)
    ax = gca ;
end
if ~exist('conn_clust_struct','var') || isempty(conn_clust_struct)
    % if we don't get cluster struct as input, load default
    rootPath = 'D:\Fly Imaging\Erica Interneuron Stacks\connectivity\' ;
    clustFn = 'conn_clust_struct.mat' ;
    
    % load clustering struct
    clustPath = fullfile(rootPath, clustFn) ;
    conn_clust_struct = importdata(clustPath) ;
end
if ~exist('clusterNum','var') || isempty(clusterNum)
    % just make "1" the default cluster number
    clusterNum = 1 ;
end
if ~exist('useSubTypeFlag','var') || isempty(useSubTypeFlag)
    % should we care about hemilineage subtype?
    useSubTypeFlag = false ;
end
if ~exist('titleFlag','var') || isempty(titleFlag)
    % put a title on the plot?
    titleFlag = false ;
end
if ~exist('colorMode','var') || isempty(colorMode)
    % how to color pie charts
    colorMode = 'hemilineage' ;  % 'hemilineage' | 'cluster' | 'neurotransmitter' | 'behavior'
end
if ~exist('scaleChartFlag','var') || isempty(scaleChartFlag)
    % scale area of pie chart to reflect cluster size?
    scaleChartFlag = true ;
end
if ~exist('sortWedgeFlag','var') || isempty(sortWedgeFlag)
    % put wedges of the same color next to each other?
    sortWedgeFlag = true ;
end
% ---------------
%  params
% ---------------
% plot params
cluster_cmap_name = 'Dark2' ; % cluster cmap
patch_alpha = 0.5 ; %0.3 ; % alpha value for pie chart patches
% lineStyleMain = 'none' ; % pie chart line style
% lineStyleOthers = '-' ; % pie chart line style
labelFontSize = 6 ;
titleFontSize = 8 ;
fontName = 'arial' ;

% search string to use to parse hemilineages
alpha_var = char(945) ;
beta_var = char(946) ;
gamma_var = char(947) ;
delta_var = char(948) ;
epsilon_var = char(949) ;
greek_str = strjoin({alpha_var, beta_var, gamma_var, delta_var,...
    epsilon_var},'') ;
searchExp = ['(?<NB>\d*)(?<notchType>[AB])' ...
    '(?<subType>[' greek_str ']?)' ...
    '\s*(?<neuromere>[t]\d)'] ;
%|'...
%     '(?<NB>\d*)(?<notchType>[AB])' ...
%     '\s*(?<neuromere>[t]\d)'] ;
parse_fields = {'NB', 'notchType', 'subType', 'neuromere'} ;

% ------------------------------------------------------
%% get hemilineage IDs of current cluster
% read out all cluster indices
clust_idx = [conn_clust_struct.cluster_idx] ;

% get membership count of each cluster
edges = 1:(max(clust_idx)+1) ;
clust_counts = histcounts(clust_idx,edges) ;
max_clust_count = max(clust_counts) ;

% get the cluster indices as permuted in dendrogram
[~, T_clust_perm_unique, ~] = get_cluster_outperm(conn_clust_struct) ;

% "clusterNumIdx" is the clust_idx number that corresponds to what I call
% clusterNum (since the order of the dendrogram (which I use as cluster
% number) doesn't directly correspond to clust_idx
clusterNumIdx = T_clust_perm_unique(clusterNum) ;

% find indices of neurons that match current cluster number
idx_curr = (clust_idx == clusterNumIdx) ;

% use these indices to get hemilineages
hemilineages_curr = {conn_clust_struct(idx_curr).hemilineageID} ;
neurons_curr = {conn_clust_struct(idx_curr).neuron_label} ;
IDs_curr = {conn_clust_struct(idx_curr).ID} ;
% --------------------------------------------------------
%% parse hemilineages
% initialize structure to store hemilineage info
parse_struct = struct() ;

% loop over current neurons
for k = 1:length(hemilineages_curr)
    % assign neuron name and label to struct
    parse_struct(k).neuron_label = neurons_curr{k} ;
    parse_struct(k).ID = IDs_curr{k} ;
    
    % now do hemilineage info. first, search for things that might break
    % our regexp search -- e.g. embryonic or abdominal
    if contains(hemilineages_curr{k}, 'abd','IgnoreCase',true)
        % in this case, it's abdominal
        parse_struct(k).NB = 'abd.' ;
        parse_struct(k).notchType = [] ;
        parse_struct(k).subType = [] ;
        parse_struct(k).neuromere = [] ;
    elseif contains(hemilineages_curr{k}, 'emb','IgnoreCase',true)
        % in this case, it's embryonic
        parse_struct(k).NB = 'emb.' ;
        parse_struct(k).notchType = [] ;
        parse_struct(k).subType = [] ;
        parse_struct(k).neuromere = [] ;
    else
        % if it's not one of these two edge cases, use regexp to find
        % detailed hemilineage info
        tokenNames = regexp(hemilineages_curr{k}, searchExp,'names') ;
        for m = 1:length(parse_fields)
            parse_struct(k).(parse_fields{m}) = tokenNames.(parse_fields{m}) ;
        end
    end
end

% --------------------------------------------------------
%% process hemilineage info for plot
% get all hemilineages from parse_struct -- include subtype?
if useSubTypeFlag
    hl_for_plot = arrayfun(@(x) strcat(x.NB, x.notchType, x.subType), ...
        parse_struct, 'UniformOutput', 0) ;
else
    hl_for_plot = arrayfun(@(x) strcat(x.NB, x.notchType), ...
        parse_struct, 'UniformOutput', 0) ;
end


% now get the unique hemilineages in current cluster, as well as their
% count
hl_unique = unique(hl_for_plot,'stable') ;
hl_count = cellfun(@(y) sum(ismember(hl_for_plot, y)), hl_unique) ;

% ---------------------------------------------------------
%% color map for plot
% get main cluster color
cluster_colors = brewermap(length(T_clust_perm_unique), ...
    cluster_cmap_name) ;
clusterColor = cluster_colors(clusterNum,:) ;

% figure out slice color depending on which mode we've selected
switch colorMode
    case 'cluster'
        % other hemilineages get lighter versions of color
        cmap_mat = [clusterColor ; 1, 1, 1 ] ;
        otherColors = customcolormap([0, 1], cmap_mat, length(hl_unique)+2) ;
        otherColors = otherColors(1:end-3,:) ;
        
        % combine into one color matrix
        pieColorMat = [clusterColor ; otherColors] ;
        pieColorMat = circshift(pieColorMat, (max_count_ind - 1), 1) ;
        
    case {'hemilineage', 'behavior', 'neurotransmitter'}
        % get current hemilineages (we should have them already, but
        % subtypes could mess this up)
        hl_for_color = arrayfun(@(x) strcat(x.NB, x.notchType), ...
            parse_struct, 'UniformOutput', 0) ;
        hl_color_unique = unique(hl_for_color,'stable') ;
        
        % get colors corresponding to hemilineage
        [pieColorMat, colorKey] = getHemilineageColor(hl_color_unique, ...
            colorMode) ;
        
    otherwise
        fprintf('Error: invalid color mode selection: %s \n', colorMode)
        keyboard
end

% ----------------------------------------------------------
%% make pie chart
% labels for pie chart -- this should be hemilineage + fraction
% pie_labels = cellfun(@(x,y) sprintf('%s\n(%d/%d)',x, y, sum(hl_count)), ...
%     hl_unique, num2cell(hl_count), 'UniformOutput', false) ;
pie_labels = cellfun(@(x,y) sprintf('%s\n(%d/%d)',x, y, sum(hl_count)), ...
    hl_unique, num2cell(hl_count), 'UniformOutput', false) ;

% sort wedges so same colors are grouped together
if sortWedgeFlag && exist('colorKey','var')
    % sort color keys
    [~, sort_idx] = sortrows(colorKey, 'ascend') ; 
    
    % apply same sorting to colormat, names, and counts
    pieColorMat = pieColorMat(sort_idx,:) ; 
    hl_count = hl_count(sort_idx) ; 
    pie_labels = pie_labels(sort_idx) ; 
    
end

% get index for maximum hemilineage count -- we'll bold the label 
[max_hl_count, max_count_ind] = max(hl_count) ;

% draw pie chart
p = pie(ax, hl_count, pie_labels) ;

% get pie chart patches and text
p_patches = findobj(p, 'Type','patch') ;
p_labels = findobj(p, 'Type','Text') ;

% scale pie chart area by cluster size?
if scaleChartFlag
    scaleFactor = sum(idx_curr)/max_clust_count ;
else
    scaleFactor = 1.0 ;
end

% loop through slices and adjust coloring/alpha/linestyle
for q = 1:length(p_patches)
    p_patches(q).FaceColor = pieColorMat(q,:) ;
    p_patches(q).FaceAlpha = patch_alpha ;
    %     if q == max_count_ind
    %         p_patches(q).LineStyle = lineStyleOthers ;
    %     else
    %         p_patches(q).LineStyle = lineStyleOthers ;
    %     end
    
    % do scaling (will have no effect if ~scaleChartFlag)
    p_patches(q).XData =  sqrt(scaleFactor).*p_patches(q).XData ;
    p_patches(q).YData =  sqrt(scaleFactor).*p_patches(q).YData ;
    p_patches(q).ZData =  sqrt(scaleFactor).*p_patches(q).ZData ;
end
% loop through labels to change font properties
for q = 1:length(p_labels)
    p_labels(q).FontName = fontName ;
    p_labels(q).FontSize = labelFontSize ;
    
    % label color also depends on colorMode
    switch colorMode
        case 'cluster'
            p_labels(q).Color = clusterColor ;
        case {'hemilineage', 'behavior', 'neurotransmitter'}
            p_labels(q).Color = pieColorMat(q,:) ;
    end
    p_labels(q).HorizontalAlignment = 'center' ;
    
    if (q == max_count_ind) && (max_hl_count > 1)
        p_labels(q).FontWeight = 'bold' ;
    end
    
    % make sure labels are on outside of circle
    label_pos = p_labels(q).Position ; % positon of text. center maybe?
    
    % do scaling if we're adjusting area
    label_pos = sqrt(scaleFactor).*label_pos ; 
    p_labels(q).Position = label_pos ; 
    
    % get extent of label
    label_extent = p_labels(q).Extent ; % extent of rectangle surrounding text
    
    % use extent to get corners -- make sure they're all outside circle
    label_corners = [label_extent(1), label_extent(2) ; ...
        label_extent(1)+label_extent(3), label_extent(2) ; ...
        label_extent(1), label_extent(2)+label_extent(4) ; ...
        label_extent(1)+label_extent(3), label_extent(2)+label_extent(4)] ;
    label_corner_dist = sqrt(dot(label_corners, label_corners,2)) ;
    
    % if any are inside circle (radius "pie_rad"), move them
    if any(label_corner_dist < sqrt(scaleFactor))
        % get corner that's closest to the center (nned to move it outside
        % circle)
        [min_corner_dist, min_ind] = min(label_corner_dist) ;
        
        % normalized position of innermost corner
        label_corner_norm = label_corners(min_ind,:) ./ ...
            norm(label_corners(min_ind,:)) ;
        
        move_amt = 1.1*sqrt(scaleFactor).*label_corner_norm - ...
            label_corners(min_ind,:) ; 
%         move_amt = 1.1*label_corner_norm - ...
%             label_corners(min_ind,:) ;
        %label_dist = norm(label_pos) ;
        p_labels(q).Position(1:2) = label_pos(1:2) + move_amt ;
    end
end

% -----------------------------------------------
% axis properties
axis(ax, 'equal') % not sure if this is necessary
if titleFlag
    h_title = title(ax, ['Cluster ' num2str(clusterNum)],...
        'FontName',fontName, ...
        'FontSize', titleFontSize,...
        'FontWeight', 'normal') ;
    h_title.Color = 'k' ; % clusterColor ;
    h_title.Position(1) = h_title.Position(1) - 1 ;
end
end