% -------------------------------------------------------------------------
% script to take a given cluster of vnc neurons, and plot them on a
% template vnc to demonstrate what they look like
% -------------------------------------------------------------------------
rootPath = 'D:\Fly Imaging\Erica Interneuron Stacks\' ;
vncPath = fullfile(rootPath, 'VNC_for_drawings') ; % need to get this!
clustPath = fullfile(rootPath, 'connectivity') ;
savePath = ['D:\Dropbox\Paper Manuscripts\Janelia vnc paper\' ...
    'figure drafts\clustering_figure\vnc_drawings'] ;
dataNameListIN = {'IN'}; % {'IN','VUM'};

% save figure?
saveFlag = true ;
saveMultiViewFlag = true ;

% skeletonize neuron images?
skelFlag = false ;
symFlag = true ;

% draw vnc on background?
vncFlag = true ;

% which cluster to look at
clusterNum = 2 ;

% plot params
% figPosition = [1921, 41, 1920/5, 963] ;
cmap_name = 'Dark2' ;
colorMode = 'cmap_same' ; % 'cmap_all' = all neurons get colored a shade of
%      the main color
% 'cmap_one' = most neurons colored gray; one is
%      given the main color
% 'cmap_same' = all neurons get main color
xlim = [90, 490] ;
ylim = [1, 1030] ;
zlim = [20, 315] ;

lineWidth = 0.5 ;
lineAlpha = 0.75 ;  % 0.5
smoothSkelFlag = false ;

vnc_color = 0.7*[1, 1, 1] ;
vnc_alpha = 0.05 ;
% -------------------------------------------------------------------
%% load data
% clustering results
clustFilename = 'conn_clust_struct.mat' ;
conn_clust_struct = importdata(fullfile(clustPath, clustFilename)) ;

% directories with neuron names
dataPathListIN = cell(length(dataNameListIN),1) ;
for k = 1:length(dataNameListIN)
    dataNameIN = dataNameListIN{k} ;
    
    % data paths depend on what we'd like to plot (symmetrized? skeleton?)
    if symFlag && skelFlag
        folderNameCurr = 'skeletonized_sym' ;
    elseif ~symFlag && skelFlag
        folderNameCurr = 'skeletonized_non_sym' ;
    elseif symFlag && ~skelFlag
        folderNameCurr = 'binarized_new' ;
    else
        folderNameCurr = 'binarized_non_sym' ;
    end
    
    % add full data path to cell array
    dataPathListIN{k} = fullfile(rootPath, dataNameIN, folderNameCurr) ;
end

% get data file suffix to look for
if skelFlag
    searchSuffix = '*skel_struct.mat' ;
else
%     searchSuffix = '*bw.mat' ;
    searchSuffix = '*bw_coords.mat' ;
end

% using path and file search suffix from above, get directory for data
% from all interneurons (includes "IN" and "VUM")
for m = 1:length(dataNameListIN)
    dataPathIN = dataPathListIN{m} ;
    dataDirIN_curr = dir(fullfile(dataPathIN, searchSuffix )) ;
    if (m == 1)
        dataDirIN = dataDirIN_curr ;
    else
        dataDirIN = vertcat(dataDirIN, dataDirIN_curr) ;
    end
end
N_IN = length(dataDirIN) ;

% *** VNC PATH INFO HERE
if vncFlag
    vnc_struct = importdata(fullfile(vncPath,'vnc_struct.mat')) ;
    vncBW = vnc_struct.vncBW ;
else
    vncBW = [] ;
end

% -------------------------------------------------------------------
%% get IN and VUM labels
label_cell_IN = cell(N_IN,1) ;
for k = 1:N_IN
    fn = dataDirIN(k).name ;
    fn_split = strsplit(fn,'_') ;
    label_str = fn_split{1} ; % strjoin(fn_split(1:end-2),'_') ;
    
    label_cell_IN{k} = label_str ;
end

% ----------------------------------------------------------------------
%% find IN (or VUM) stacks that belong to current cluster
% (also get color for current cluster)
clust_idx = [conn_clust_struct.cluster_idx] ;
clust_idx_unique = unique(clust_idx) ;

cluster_colors = brewermap(length(clust_idx_unique),cmap_name) ;

% indices for cluster from clustering struct -- note that cluster number in
% conn_clust_struct does not refer to cluster numbers as I've defined them
% (order in the dendrogram plot). so need to define clusterNumIdx to match
% them up
[~, T_clust_perm_unique, ~] = get_cluster_outperm(conn_clust_struct) ;
clusterNumIdx = T_clust_perm_unique(clusterNum) ;

% find indices of neurons that match current cluster number
idx_curr = find(clust_idx == clusterNumIdx) ;

% color map for plot
clusterColor = cluster_colors(clusterNum,:) ;

% corresponding indices for directory
idx_for_dir = zeros(size(idx_curr)) ;
for q = 1:length(idx_curr)
    idxx = idx_curr(q) ;
    annotation_label = conn_clust_struct(idxx).neuron_label ;
    annotation_label = strrep(annotation_label,'_','-') ;
    mat_match_idx = cellfun(@(y) strcmp(y, annotation_label),...
        label_cell_IN) ;
    try
        idx_for_dir(q) = find(mat_match_idx) ;
    catch
        keyboard
    end
end

dataDirIN = dataDirIN(idx_for_dir) ;
N_IN = length(dataDirIN) ;

% -------------------------------------------------------
%% define custom colormap based on current cluster color
% random number generator seed (for color shuffling)
rng(47)
switch colorMode
    case 'cmap_all'
        % generate colormap based on primary figure color
        cmap_mat = [1,1,1 ; clusterColor ; 0,0,0 ] ;
        clusterColorCmap = customcolormap([0, 0.75, 1], cmap_mat, ...
            N_IN+4) ;
        neuronColors = clusterColorCmap(3:end-1,:) ;
        
        % then shuffle colors so we don't get things right next to each
        % other with same color
        
        perm_ind = randperm(size(neuronColors,1)) ;
        neuronColors = neuronColors(perm_ind,:) ;
    case 'cmap_one'
        perm_ind = randperm(N_IN) ;
        neuronColors = zeros(N_IN,3) ;
        neuronColors(perm_ind(1:2),:) = repmat(clusterColor,2,1) ;
        neuronColors(perm_ind(3:end),:) = ...
            repmat(linspace(0.6,0.9,N_IN-2)',1,3).*ones(N_IN-2,3) ;
    case 'cmap_none'
        neuronColors = (0.8).*ones(N_IN,3) ;
    case 'new_cmap'
        neuronColors = viridis(N_IN) ; % brewermap(N_IN,'Spectral') ;
        perm_ind = randperm(size(neuronColors,1)) ;
        neuronColors = neuronColors(perm_ind,:) ;
    case 'cmap_same'
        neuronColors = repmat(clusterColor, N_IN, 1) ;
end

% ----------------------------------------------------------------------
%% loop through and load neuron data to plot
% initialize figure and axis
h_main = figure ;
ax = gca ;
hold on

% print update
fprintf('Drawing neurons...\n')

% loop over IN data structures and draw/save
for k = 1:N_IN
    
    % create a parent hgtform to hopefully use to move stuff around later
    parent = hgtransform('Parent',ax);

    % plot vnc surface (if plotting skeleton; otherwise the patch code draws it
    % already)
    %     if vncFlag && skelFlag
    %         ax = drawVNCBackgroundToAxis(ax, h_main, vnc_struct, vnc_color,...
    %             vnc_alpha ) ;
    %     end
    if vncFlag
        ax = drawVNCBackgroundToAxis(ax, h_main, vnc_struct, vnc_color, ...
            vnc_alpha, parent) ;
    end
    
    
    % get IN identity
    fn = dataDirIN(k).name ;
    [di1, di2] = regexp(fn, 'SS\d{5}') ;
    driver = fn(di1:di2) ;
    
    fn_split = strsplit(fn,'_') ;
    line_name = fn_split{1} ;
    
    % load interneuron
    IN_path_curr = fullfile(dataDirIN(k).folder, fn) ;
    IN_BW = importdata(IN_path_curr) ;
    
    % draw neuron and vnc(?)
    if ~skelFlag
        % draw neuron
        ax = drawNeuronToAxis(ax, IN_BW, neuronColors(1,:), [], [], parent) ;
        % ax = drawNeuronsAndVNC(ax, h_main, vncBW, bwCell, neuronColors) ;
    else
        for q = 1:length(IN_BW)
            ax = drawSkelGraph(ax, IN_BW(q), neuronColors(k,:), lineWidth, ...
                lineAlpha, smoothSkelFlag, ~vncFlag) ;
        end
    end
    
    % ------------------------------
    % set axis properties
    set(ax,'DataAspectRatio',[1 1 1])
    set(ax,'PlotBoxAspectRatio',[1 1 1])
    
    % ----------------------------------------
    % save results?
    if saveFlag
        fprintf('Saving neuron image %d/%d... \n', k, N_IN)   
        
        if skelFlag
            set(h_main, 'Renderer','painters')
            savePrefix = 'skel' ;
            % also need to rotate view if this is the case
            view(180, 90) ; % this gets it to horizontal view
        else
%             % add lighting?
%             camlight;
%             lighting phong
            savePrefix = 'bw' ;
        end
        
        % either save three views of VNC or just horizontal view
        if saveMultiViewFlag
            % make sure camera doesn't change zoom
            camva('manual')
            
            % list of 3 primary views for vnc
            view_list = {[180, 90], [0, 0], [90, 0]} ; % horizontal, transverse (coronal), sagittal
            view_names = {'horz', 'trans', 'sag'} ;
            
            % loop over views and save plots
            for kk = 1:length(view_list)
                % set current view
                view_curr = view_list{kk} ;
                view(view_curr(1), view_curr(2)) ;
                
                % adjust lighting for new view (if not plotting skels)
                if ~skelFlag
                    camlight;
                    lighting phong
                end
                
                % generate save filename
                saveFn = strjoin({savePrefix, 'neurons', 'cluster', ...
                    num2str(clusterNum,'%02d'), 'neuron', num2str(k,'%02d'),...
                    line_name, view_names{kk}}, '_') ;
                savePathCurr = fullfile(savePath, saveFn) ;
                
                % save png of image
                export_fig(h_main, savePathCurr, '-dpng', '-r500')
            end
        else
            % if only saving one view, do that here
            saveFn = strjoin({savePrefix, 'neurons', 'cluster', ...
                num2str(clusterNum,'%02d'), 'neuron', num2str(k,'%02d'),...
                line_name}, '_') ;
            savePathCurr = fullfile(savePath, saveFn) ;
            export_fig(h_main, savePathCurr, '-dpng', '-r500')
        end
    end
    
    % print updates
    fprintf('Completed %d/%d lines \n',k,N_IN)
    
    % clear axis
    cla(ax)
end


