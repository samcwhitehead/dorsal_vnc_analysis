% -------------------------------------------------------------------------
% script to generate multipe matrix + dendro plots for annotated anatomy
% -------------------------------------------------------------------------
% path info 
[mfilePath, ~, ~] = fileparts(mfilename('fullpath')) ; 
figDirectory = fileparts(mfilePath) ;
parentDirectory = fileparts(figDirectory) ; 
rootPath = fullfile(parentDirectory, 'data', 'clustering') ;
savePath = fullfile(mfilePath,'output') ;

% ----------------------------------
% options and params
saveFlag = false ; % for plots
colorClustFlag = false ;
excludeVUMFlag = true ;

plotMode = 'all' ; % 'all' | 'individual' (give each cluster its own figure or plot all together)
colorPatchFlag = true ; % color patches for cluster boundaries?
flipYTickLabelFlag = true ; % flip labels, so that when we rotate a horizontal plot they look right

DNFlag = false ; % make plot for DNs or INs + VUMs (default)?
newAnnotFlag = true ; % use new, more precise annotations?

capitalScale = 2 ;  % relative weighting of strong vs weak i/o
useIpsiContraFlag = true;  % distinguish between ipsi and contra sides? 
rmWeakAnnotFlag = true;  % remove weak i/o annotations?
onlyMatFlag = true ; % stop after plotting matrix because digraph code still needs updating
colorIpsiContraFlag = true ;  % color ipsi- and contralateral labels differently?

% switch filename for data to load depending on neuron type
if DNFlag
    suffixStr = 'DN' ;
    fprintf('Not working at the moment \n')
    keyboard
else
    suffixStr = 'IN' ;
end
dataFilename = sprintf('conn_clust_struct_%s.mat', suffixStr) ;

% make sure save directory exists
if ~exist(savePath,'dir')
    mkdir(savePath)
end

% if using new annotations and distinguishing between ipsi and
% contralateral annotations
side_list = {'ipsi', 'contra'} ;

% playing around with coloring ipsi and contrlateral labels differently
ipsiColor = [0.0, 0.0, 0.0] ;
contraColor = 0.5*[1.0, 1.0, 1.0] ;

% --------------------------------------------------
% digraph params
layout = 'circle' ; % 'custom' ; % % was circle %'custom' | 'circle' | 'custom circle' | 'circle and patch'
binEdgeFlag = true ;
if ismember(layout, {'force','circle and patch'})
    rmNodesFlag = true ;
else
    rmNodesFlag = false ; % false
end

normEdgeFlag = false ; % normalize edge widths to current max?

%------------------------------------
% plot params
figUnits = 'inches' ;
switch plotMode
    case 'all'
        if DNFlag
            figPositionClust = [0.5, 3.4896, 5.826, 3.7] ; 
            figPositionGraph = [0.5, 2.3229, 4.9, 5.6667] ;  
        else
            figPositionClust = [0.5, 3.4896, 8.8, 3.2] ; 
            figPositionGraph = [0.5000, 0.7083, 5.4, 7.2813] ;
            if useIpsiContraFlag
                figPositionClust(4) = 2.4*figPositionClust(4) ;
            end
            
        end
    case 'individual'
        figPositionClust = [0.5, 3.4896, 19.6667, 4.3958] ;
        if useIpsiContraFlag
            figPositionGraph = [0.7500, 4.1250, 3.2500, 3.5417] ;
        else
            figPositionGraph = [0.5, 4.1250, 6.8, 2.1042] ;  
        end
        
    otherwise
        fprintf('Invalid plot mode: %s \n',plotMode)
        keyboard
end
%figPosition = [20.2188    3.4896   19.6667    3.3646] ;
axisFontSize = 6 ;
axisFontSizeSmall = 4.5 ;

% axis position (dendrogram and matrix):
gap = [0.01, 0.01] ; %[0.1 0.10] ;
marg_h = [0.26, 0.26] ; %[0.1 0.02] ;
if DNFlag
    marg_w = [0.075, 0.001] ;
else
    marg_w = [0.04 0.001] ;
end

n_rows = 5 ; % 4
n_cols = 1 ;

% axis position (digraphs):
gap_g = [0.04, 0.075] ; %[0.03, 0.05] ; %[0.1 0.10] ;
marg_h_g = [0.03, 0.03] ; %[0.08 0.08] ;
marg_w_g = [0.07, 0.07] ; %[0.035, 0.03] ;

n_rows_g = 5 ;
n_cols_g = 4 ;

xlim_g = 1.75*[-1,1] ; % [-1.5   1.5] ;
ylim_g = 1.75*[-1,1] ; % [-1.5   1.5] ;

xlim_vnc = [-91, 1600] ;
ylim_vnc = [1, 425] ;

% ------------------------------------------------------------------------
%% load previous clustering results
connPath = fullfile(rootPath, dataFilename) ;
if exist(connPath,'file')
    conn_clust_struct = importdata(connPath) ;
else
    fprintf('No clustering results found -- quitting \n')
    return
end

% read in clustering and annotation results
T_clust = [conn_clust_struct.cluster_idx]' ;
tree  = conn_clust_struct(1).tree ;
neuron_labels = {conn_clust_struct.neuron_label} ;

% read ou annotation values. each neuropil target has its own struct field,
% so need to find and read out each
conn_clust_fields = fieldnames(conn_clust_struct) ;
target_idx = cellfun(@(y) contains(y, 'target_','IgnoreCase',true), ...
    conn_clust_fields) ;
target_ind = find(target_idx) ;

% initialize output
if newAnnotFlag && useIpsiContraFlag
    % check whether or not we're using ipsi and contra separately
    annotationVals = nan(length(neuron_labels), 2*length(target_ind)) ;
    target_labels = cell(2*length(target_ind),1) ;
    
else
    % otherwise just use normal
    annotationVals = nan(length(neuron_labels), length(target_ind)) ;
    target_labels = cell(length(target_ind),1) ;
    
end


% loop through target labels
for k = 1:length(target_ind)
    % get field of conn_clust_struct corresponding to current target
    % neuropil
    ind_curr = target_ind(k) ;
    field_curr = conn_clust_fields{ind_curr} ;
    
    % get label of neuropil from field name
    targ_label = regexp(field_curr,'[^target_]\w*','match') ;
    
    % read out current annotation values
    annotValsCurr = [conn_clust_struct.(field_curr)] ;
    
    % add label and annotation vals to arrays initialized above
    if useIpsiContraFlag && newAnnotFlag
        % if separating ipsi and contra, add these as separate entries
        ipsi_ind = 2*(k-1) + 1 ;
        contra_ind = 2*k ;
        
        % labels
        target_labels{ipsi_ind} = sprintf('%s_{il}',targ_label{1}) ;
        target_labels{contra_ind} = sprintf('%s_{cl}',targ_label{1}) ;
        
        % annotation values -- need to loop over to convert the 1x4 array
        % to something in the 1-3ish range. this is going to be annoying.
        % to mitigate annoyance
        for m = 1:size(annotValsCurr,2)
            % current neuron annotations
            neuronAnnot = annotValsCurr(:,m) ;
            % ipsilateral value
            annotationVals(m,ipsi_ind) = ...
                newAnnotationVecToNum(neuronAnnot([1,3]), capitalScale,...
                [], rmWeakAnnotFlag) ;
            % contralateral value
            annotationVals(m,contra_ind) = ...
                newAnnotationVecToNum(neuronAnnot([2,4]), capitalScale,...
                [], rmWeakAnnotFlag) ;
        end
        
    elseif size(annotValsCurr,1) > 1
        % in this case, we're using new annotations, but need to collapse
        % into one number per (neuron, neuropil) combo
        target_labels{k} = targ_label{1} ;
        
        % do the same loop over neurons as above to collapse down to 1,2,3
        for m = 1:size(annotValsCurr,2)
            neuronAnnot = annotValsCurr(:,m) ;
            annotationVals(m,k) = ...
                newAnnotationVecToNum(neuronAnnot, capitalScale,...
                [], rmWeakAnnotFlag) ;
        end
        
    else
        % otherwise just read directly in
        target_labels{k} = targ_label{1} ;
        annotationVals(:,k) = annotValsCurr;
    end
end

% resort ipsi and contralateral so they form two separate column blocks
if useIpsiContraFlag
    resort_ind = [1:2:length(target_labels), 2:2:length(target_labels)] ;
    annotationVals = annotationVals(:, resort_ind) ;
    target_labels = target_labels(resort_ind) ;
end

% -------------------------------------------------------------------
%% need to read out some other things if we're doing individual plots
if strcmp(plotMode, 'individual')
    % clustering params
    linkageMethod = conn_clust_struct(1).linkageMethod ;
    distanceType = conn_clust_struct(1).distanceType ;
    
    % also need reshaped annotation matrix
    if newAnnotFlag
        % if we're using new annotations, easiest for now to just re-do
        % conversion.
        % NB: for ease of graph plotting, we probably want to set
        % capitalScale to 1 and collapse ipsi/contralateral onto each other
        annotationValsCell = conn_clust_struct(1).annotationValsCell ;
        annotationValsReshape = ...
            convertNewAnnotationVals(annotationValsCell,1,false) ;
        
    else
        annotationValsReshape = reshapeAnnotationMat(annotationVals) ;
    end
end
% -------------------------------------------------------------------------
%% get things ready to plot
% flip matrices and labels
annotationValsFlip = fliplr(annotationVals) ;
target_labels_flip = fliplr(target_labels) ;

% color info
colorThresh = nan ; % tree(end-max_clusts+2,3)-eps + 0.002; % color threshold for dendrogram
if DNFlag
    % colormap  for *DN* clusters
    %     cluster_cmap = 'Set1' ;
    pt_blue = [0, 68, 136]./255 ;
    pt_red = [187, 85, 102]./255 ;
    pt_yellow = [221, 170, 51]./255 ;
    
    pt_darkblue = [34, 34, 85]./255 ;
    pt_darkcyan = [34, 85, 85]./255 ;
    pt_darkgreen = [34, 85, 34]./255 ;
    
    pt_indigo = [51, 34, 136]./255 ;
    pt_teal = [68, 170, 153]./255 ;
    pt_olive = [153, 153, 51]./255 ;
    pt_wine = [136, 34, 85]./255 ;
    
    % grays
    light_gray = 0.8*[1, 1, 1] ;
    medium_gray = 0.6*[1, 1, 1] ;
    dark_gray = 0.4*[1, 1, 1] ;
    
    %     cluster_colors = [pt_yellow ; pt_red ; pt_blue] ;
    %     cluster_colors = [pt_indigo ; pt_olive ; pt_wine] ;
    cluster_colors = [medium_gray ; dark_gray ; light_gray] ;
    cluster_cmap = cluster_colors ;
else
    % colormap *IN* for clusters
    cluster_cmap = 'Dark2' ;
    cluster_colors = brewermap(50,cluster_cmap) ;
end


% unique cluster indices
T_unique = unique(T_clust) ;

%  get axis limits for full data
xlim = [0.5, length(neuron_labels) + 0.5] ;
ylim2 = [0.5, length(target_labels) + 0.5] ;
ylim1 = [0, max(tree(:,3))] ;
ylim_list = {ylim1, ylim2} ;


% ---------------------------------------------------------------
%% make figure window and subplot axes for matrix + dendrogram
% IF PLOTTING ALL
if strcmp(plotMode,'all')
    h_clust = figure('PaperPositionMode','auto','MenuBar','none',...
        'ToolBar','none','DockControls','off','Units',figUnits,...
        'OuterPosition',figPositionClust) ;
    
    % initialize axes (ax1 for dendrogram, ax2 for matrix)
    ax1 = subtightplot(n_rows, n_cols, 1, gap, marg_h, marg_w, ...
        'Parent', h_clust) ;
    ax2 = subtightplot(n_rows, n_cols, 2:n_rows, gap, marg_h, marg_w, ...
        'Parent',h_clust) ; % 2:n_rows
    
    % populate axes
    [ax1, ax2, outperm] = plot_connectivity_mat(h_clust, ax1, ax2,  ...
        annotationValsFlip, tree,  neuron_labels, ...
        target_labels_flip, T_clust, colorThresh, [0, 0, 0], ...
        colorClustFlag, xlim,  ylim_list, axisFontSizeSmall,...
        [], colorPatchFlag, flipYTickLabelFlag, cluster_cmap) ;
    
    % color ipsi and contralateral tick labels differently?
    if useIpsiContraFlag && colorIpsiContraFlag
        % read out current tick labels and make a copy for colored version
        np_tick_labels = get(ax2, 'YTickLabel') ;
        np_tick_labels_new = np_tick_labels ;
        
        % loop over tick labels and change their color depending on
        % subscript
        for q = 1:length(np_tick_labels)
            tick_label = np_tick_labels{q} ;
            
            % change color depending on subscript
            if contains(tick_label, '_{il}')
                np_tick_labels_new{q} = sprintf('{{\\color[rgb]{%f, %f, %f} %s}}',...
                    ipsiColor(1), ipsiColor(2), ipsiColor(3), tick_label) ;
            elseif contains(tick_label, '_{cl}')
                np_tick_labels_new{q} = sprintf('{{\\color[rgb]{%f, %f, %f} %s}}',...
                    contraColor(1), contraColor(2), contraColor(3), ...
                    tick_label) ;
            end
        end
        
        % assign new tick labels
        set(ax2, 'YTickLabel', np_tick_labels_new) ;
    end
    % save results?
    if saveFlag
        % matrix
        matSaveFn = sprintf('connectivity_cluster_all_%s',suffixStr) ;
        %         export_fig(h_clust,fullfile(savePath, ...
        %             [matSaveFn '.png']),'-dpng','-r600')
        exportgraphics(h_clust,fullfile(savePath, [matSaveFn '.png']),...
            'Resolution',600)
        print(h_clust,fullfile(savePath, [matSaveFn '.svg']),'-dsvg')
    end
end

% if we're only clustering DNs, stop here (no inputs, so no need for graph
if (DNFlag || onlyMatFlag) && ~strcmp(plotMode,'individual')
    return
end
% ---------------------------------------------------------------
%% make figure window and subplot axes for digraph
N_clusts = length(T_unique) ;
[T_clust_perm, T_clust_perm_unique, outperm] = ...
    get_cluster_outperm(conn_clust_struct) ;
% T_clust_perm = T_clust(outperm) ;
% T_clust_perm_unique = unique(T_clust_perm, 'stable') ;

% ax_array = gobjects(N_clusts,1) ;
if strcmp(plotMode,'all')
    h_graph = figure('PaperPositionMode','auto','MenuBar','none',...
        'ToolBar','none','DockControls','off','Units',figUnits,...
        'OuterPosition',figPositionGraph) ;
end

% ---------------------------------------------------------------------
% before making plots for individual clusters, modify target labels for
% plots:
target_labels_mod = target_labels ;

% loop over target labels
for m = 1:length(target_labels_mod)
    if ~newAnnotFlag && contains(target_labels_mod{m},...
            'LegNp','IgnoreCase',1)
        % shorten some target label abbreviations for things to fit
        target_labels_mod{m} = target_labels_mod{m}(1:end-2) ;
    end
    %      elseif newAnnotFlag && useIpsiContraFlag
    %          % remove underscores since, for graph plots, we need to collapse
    %          % ipsi/contra distinction
    %          target_labels_mod{m} = target_labels_mod{m}(1:end-5) ;
    %      end
    if contains(target_labels_mod{m}, '_{il}')
        target_labels_mod{m} = sprintf('{{\\color[rgb]{%f, %f, %f} %s}}',...
            ipsiColor(1), ipsiColor(2), ipsiColor(3), target_labels_mod{m}) ;
    elseif contains(target_labels_mod{m}, '_{cl}')
        target_labels_mod{m} = sprintf('{{\\color[rgb]{%f, %f, %f} %s}}',...
            contraColor(1), contraColor(2), contraColor(3), ...
            target_labels_mod{m}) ;
    end
end
target_labels_mod = unique(target_labels_mod, 'stable') ;

% ----------------------------------------------------------------------
% get annotation values specifically to be used for directed graphs
annotationValsGraph = round(annotationVals) ;
target_labels_graph = target_labels_mod ;

% flip contra order to make symmetric
if useIpsiContraFlag
    % get flip index order
    label_sort_ind = 1:length(target_labels_graph) ;
    contra_ind = (length(target_labels_graph)/2 + 1): ...
        length(target_labels_graph) ;
    label_sort_ind(contra_ind) = fliplr(label_sort_ind(contra_ind)) ;
    
    % apply new sorting
    annotationValsGraph = annotationValsGraph(:, label_sort_ind) ;
    target_labels_graph = target_labels_graph(label_sort_ind) ;
end

% remove nodes that never get used
no_annot_idx = (sum(annotationValsGraph,1) < 1) ;
% ... but try to keep some symmetry is using ipsi/contra
if useIpsiContraFlag
    N_labels = length(target_labels_graph) ;
    no_annot_idx_sym = no_annot_idx(1:(N_labels/2)) & ...
        fliplr(no_annot_idx((N_labels/2 + 1):end)) ;
    no_annot_idx = [no_annot_idx_sym, fliplr(no_annot_idx_sym)] ;
end
annotationValsGraph = annotationValsGraph(:,~no_annot_idx) ;
target_labels_graph = target_labels_graph(~no_annot_idx) ;

% shift graph labels/annotations so that ipsi- and contralateral sides go
% on left and right of circle
idx_shift = 6 ; 
target_labels_graph = circshift(target_labels_graph, idx_shift) ; 
annotationValsGraph = circshift(annotationValsGraph, idx_shift, 2) ; 
% --------------------------------------------
% loop over clusters
for ind = 1:N_clusts
    % find indices for neurons belonging to current cluster
    idx_curr = (T_clust == T_clust_perm_unique(ind)) ;
%     if sum(idx_curr) < 2
%         keyboard
%     end
    annotationVals_curr = annotationVals(idx_curr,:) ;
    annotationValsGraph_curr = annotationValsGraph(idx_curr,:) ;
    neuron_labels_curr = neuron_labels(idx_curr) ;
    
    % get current cluster color
    colorCurr = cluster_colors(ind,:) ; % cluster_colors(T_clust_perm_unique(ind),:) ; % cluster_colors(ind,:) ;
    
    switch plotMode
        case 'all'
            % generate axis for graph plot
            ax = subtightplot(n_rows_g, n_cols_g, ind, gap_g, ...
                marg_h_g, marg_w_g) ;
            set(ax,'Parent',h_graph)
            
            % fill graph axis
            ax = plotAnnotationDigraph(ax, annotationValsGraph_curr,...
                target_labels_mod, colorCurr, layout, binEdgeFlag, ...
                [], [], rmNodesFlag, normEdgeFlag) ;
            
            % adjust axis
            if ismember(layout,{'circle','custom circle'}) && ...
                    ~(isempty(xlim_g) || isempty(ylim_g))
                set(ax,'xlim',xlim_g,'ylim',ylim_g)
            end
        case 'individual'
            % get hierarchical structure within current cluster
            if sum(idx_curr) > 1
                D_curr = pdist(annotationValsReshape(idx_curr, :), ...
                    distanceType) ;
                tree_curr = linkage(D_curr, linkageMethod) ;
                
                % ---------------------------------------------------------
                % make figure window and subplot axes for matrix and 
                % dendrogram plots. NB: only works with more than one
                % neuron
                h_clust = figure('PaperPositionMode','auto',...
                    'MenuBar','none', 'ToolBar','none',...
                    'DockControls','off','Units',figUnits,...
                    'OuterPosition',figPositionClust) ;
                
                ax1 = subtightplot(n_rows, n_cols, 1, gap, marg_h, marg_w, ...
                    'Parent', h_clust) ;
                
                ax2 = subtightplot(n_rows, n_cols, 2:n_rows, gap, marg_h,...
                    marg_w, 'Parent',h_clust) ; % 2:n_rows
                % populate figure
                plot_connectivity_mat(h_clust, ax1, ax2, ...
                    fliplr(annotationVals_curr), tree_curr, ...
                    neuron_labels_curr, fliplr(target_labels), [], ...
                    colorThresh, colorCurr, colorClustFlag, xlim, ...
                    ylim_list, axisFontSize, axisFontSize, false, false) ;
                
                % save matrix plot?
                if saveFlag
                    % matrix
                    exportgraphics(h_clust,fullfile(savePath, ...
                        ['connectivity_cluster_' num2str(ind,'%02d') ...
                        '.png']), 'Resolution', 300)
                    print(h_clust,fullfile(savePath, ...
                        ['connectivity_cluster_' num2str(ind,'%02d') ...
                        '.svg']), '-dsvg')
                end
            end
            %tree_curr = conn_clust_struct(1).tree ;

            % ---------------------------------------------------------------
            % make figure window and subplot axes for digraph
            h_graph = figure('PaperPositionMode','auto','MenuBar','none',...
                'ToolBar','none','DockControls','off','Units',figUnits,...
                'OuterPosition',figPositionGraph) ;
            ax3 = axes ; %subtightplot(n_rows, n_cols, 1:n_rows, gap, marg_h, marg_w) ;
            set(ax3,'Parent',h_graph) ;
            
            ax3 = plotAnnotationDigraph(ax3, annotationValsGraph_curr, ...
                target_labels_graph, colorCurr, layout, binEdgeFlag, ...
                [],[],rmNodesFlag, normEdgeFlag) ;
            
            % set axis limits
            if ismember(layout,{'circle','custom circle'}) && ...
                    ~(isempty(xlim_g) || isempty(ylim_g))
                set(ax3,'xlim',xlim_g,'ylim',ylim_g) ;
            elseif strcmp(layout,'custom')
                set(ax3,'xlim',xlim_vnc, 'ylim', ylim_vnc)
            end
            
            % ---------------------------------------------------------------
            % save  graph results?
            if saveFlag
                
                % digraph
                % axis(ax3,'normal') ; 
                set(h_graph, 'Renderer','painters')
                graph_fn = ['connectivity_graph_' num2str(ind,'%02d') ...
                    '_' layout] ;
                graphSavePath = fullfile(savePath, graph_fn) ;
                exportgraphics(h_graph,[ graphSavePath '.png'],...
                    'Resolution',500)
                print(h_graph,[graphSavePath '.svg'], '-dsvg')
                close all
            end
            fprintf('Completed %d/%d plots \n', ind, length(T_unique))
    end
end

% ---------------------------------------------------------------
%% save results
if saveFlag && strcmp(plotMode,'all')
    % digraph
    digraphSaveFn = sprintf('connectivity_graph_all_%s_%s',layout, suffixStr) ;
    exportgraphics(h_graph,fullfile(savePath, [digraphSaveFn '.png']),...
        'Resolution', 500)
    print(h_graph,fullfile(savePath, [digraphSaveFn '.svg']),'-dsvg')
end


