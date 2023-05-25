% -------------------------------------------------------------------------
% function to summarize hemilineage identity of interneurons by cluster in
% plot form
% -------------------------------------------------------------------------
function ax = plotClusterHemilineageBarToAxis(ax, conn_clust_struct, ...
    clusterNum, useSubTypeFlag,colorMode, sortStackFlag, ...
    fullLabelFlag, horizontalFlag)
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
if ~exist('colorMode','var') || isempty(colorMode)
    % how to color pie charts
    colorMode = 'hemilineage' ;  % 'hemilineage' | 'cluster' | 'neurotransmitter' | 'behavior'
end
if ~exist('sortStackFlag','var') || isempty(sortStackFlag)
    % put wedges of the same color next to each other?
    sortStackFlag = true ;
end
if ~exist('fullLabelFlag','var') || isempty(fullLabelFlag)
    % if true, label each stacked bar section w/ hemilineage name and
    % fraction. if false, label with just total number of cells in clust
    fullLabelFlag = true ;
end
if ~exist('horizontalFlag','var') || isempty(horizontalFlag)
    % if true, make the bars horizontal
    horizontalFlag = false ;
end
% ---------------
%  params
% ---------------
% plot params
cluster_cmap_name = 'Dark2' ; % cluster cmap
patch_alpha = 0.5 ; %0.75 ; % alpha value for pie chart patches
bar_width = 0.75 ;
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
        parse_struct(k).NB = 'abd' ;
        parse_struct(k).notchType = [] ;
        parse_struct(k).subType = [] ;
        parse_struct(k).neuromere = [] ;
    elseif contains(hemilineages_curr{k}, 'emb','IgnoreCase',true)
        % in this case, it's embryonic
        parse_struct(k).NB = 'emb' ;
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
        barColorMat = [clusterColor ; otherColors] ;
        barColorMat = circshift(barColorMat, (max_count_ind - 1), 1) ;
        
    case {'hemilineage', 'behavior', 'neurotransmitter'}
        % get current hemilineages (we should have them already, but
        % subtypes could mess this up)
        hl_for_color = arrayfun(@(x) strcat(x.NB, x.notchType), ...
            parse_struct, 'UniformOutput', 0) ;
        hl_color_unique = unique(hl_for_color,'stable') ;
        
        % get colors corresponding to hemilineage
        [barColorMat, colorKey] = getHemilineageColor(hl_color_unique, ...
            colorMode) ;
        
    otherwise
        fprintf('Error: invalid color mode selection: %s \n', colorMode)
        keyboard
end

% ----------------------------------------------------------
%% make bar for current cluster
% labels for bar plots -- this should be hemilineage + fraction

%     bar_labels = cellfun(@(x,y) sprintf('%s\n(%d/%d)',x, y, sum(hl_count)), ...
%         hl_unique, num2cell(hl_count), 'UniformOutput', false) ;
bar_labels = cellfun(@(x) sprintf('%s',x), hl_unique, ...
    'UniformOutput', false) ;


% sort wedges so same colors are grouped together
if sortStackFlag && exist('colorKey','var')
    % sort color keys (old)
    % [~, sort_idx] = sortrows(colorKey, 'ascend') ;
    
    % new attempt to try to keep bar sorting consistent across plots
    % ... first need to convert hemilineage numbers to 2 digit numbers, so
    searchStrCurr = '(?<nb_num>\d+)(?<notch_type>[AB]{1})' ; 
    hl_label_nums = zeros(size(bar_labels)) ; 
    for q = 1:length(bar_labels)
       tokenNames = regexp(bar_labels{q}, searchStrCurr, 'names') ;
       if ~isempty(tokenNames)
           hl_label_nums(q) = str2double(tokenNames.nb_num) ; 
       else
           hl_label_nums(q) = -1 ; 
       end
       
    end
    [~, sort_idx] = sort(hl_label_nums, 'ascend') ;
    
    % apply same sorting to colormat, names, and counts
    barColorMat = barColorMat(sort_idx,:) ;
    hl_count = hl_count(sort_idx) ;
    bar_labels = bar_labels(sort_idx) ;
    
end

% get index for maximum hemilineage count -- we'll bold the label
[max_hl_count, max_count_ind] = max(hl_count) ;

% normalize data for bar plot
hl_count_norm = hl_count./sum(hl_count) ;

% draw bar
if horizontalFlag
    b = barh(ax, clusterNum, hl_count_norm, bar_width, 'stacked') ;
else
    b = bar(ax, clusterNum, hl_count_norm, bar_width, 'stacked') ;
end

% get y values of bars
if length(b) == 1
    bar_y_centers = 0.5 ;
else
    bar_y_vals = [0 , cumsum([b(:).YData])];
    bar_y_centers = (bar_y_vals(1:end-1) + bar_y_vals(2:end))./2 ;
end

% loop through bar objects and adjust colors/other properties
for b_num = 1:length(b)
    % faces
    b(b_num).FaceColor = barColorMat(b_num,:) ;
    b(b_num).FaceAlpha = patch_alpha ;
    
    % edges
    b(b_num).EdgeColor = 'k' ;
    b(b_num).EdgeAlpha = patch_alpha ;
    
    % baseline
    b(b_num).BaseLine.LineStyle = 'none';
    
    % add label to each bar stack?
    if fullLabelFlag
        % bold label of hemilineage represented fraction-wise
        if (b_num == max_count_ind) && (max_hl_count > 1) 
            fontWeight = 'bold' ;
        else
            fontWeight = 'normal' ;
        end
        
        % add text
        if horizontalFlag
            text(ax, bar_y_centers(b_num), clusterNum, bar_labels{b_num}, ...
                'FontName', fontName, 'fontSize', labelFontSize, ...
                'FontWeight', fontWeight, 'Color', 'k',...
                'HorizontalAlignment','center') ;
        else
            text(ax, clusterNum, bar_y_centers(b_num), bar_labels{b_num}, ...
                'FontName', fontName, 'fontSize', labelFontSize, ...
                'FontWeight', fontWeight, 'Color', 'k',...
                'HorizontalAlignment','center') ;
        end
            
    end  
   
end

% add label at top giving full number
basic_label_str = sprintf('N=%d', sum(hl_count)) ;

% add text
if horizontalFlag
    text(ax, 1.01, clusterNum, basic_label_str, ...
        'FontName', fontName, 'fontSize', labelFontSize, ...
        'FontWeight', 'normal', 'Color', 'k',...
        'HorizontalAlignment','left',...
        'VerticalAlignment', 'middle',...
        'FontAngle','italic') ;
else
    text(ax, clusterNum, 1.02, basic_label_str, ...
        'FontName', fontName, 'fontSize', labelFontSize, ...
        'FontWeight', 'normal', 'Color', 'k',...
        'HorizontalAlignment','left',...
        'VerticalAlignment', 'bottom',...
        'Rotation',35, ...
        'FontAngle','italic') ;
end

end