% -------------------------------------------------------------------------
% function to sort neuron labels (trying to write functions to clean up
% that plotting code)
%
% INPUTS:
%   - labels: cell array of labels (strs or char arrays)
%   - neuronType: string describing the neurons we're using ('IN' | 'DN' |
%       'VUM' | 'MN' )
%   - useINClusterFlag: boolean. sort using connectivity clusters?
%   - MN_order: cell array giving motor neuron order
%   - VUM_order: cell array giving VUM neuron order
%   - conn_clust_struct: structure containing info on connectivity
%       clustering
% -------------------------------------------------------------------------
function [labels_sort, sort_ind] = sort_neuron_labels_overlap(labels, ...
    neuronType, useINClusterFlag, useDNClusterFlag, MN_order, VUM_order,...
    conn_clust_struct)
% ------------------------------------
%% inputs and params
if ~exist('useINClusterFlag','var') || isempty(useINClusterFlag)
    useINClusterFlag = true ;
end
if ~exist('useDNClusterFlag','var') || isempty(useDNClusterFlag)
    useDNClusterFlag = true ;
end
if ~exist('MN_order','var') || isempty(MN_order)
    % default order for plotting MN data
    MN_order = {'DLM', 'DVM', 'tpN', 'tp1', 'tp2', 'ps1', 'i1', 'i2', ...
        'hg1', 'hg2'} ;
end
if ~exist('VUM_order','var') || isempty(VUM_order)
    % default order for plotting MN data
    VUM_order = {'V01', 'V12', 'V03', 'V31', 'V33', 'mVUM201', 'mVUM202', ...
        'mVUM203', 'mVUM204', 'mVUM205', 'V18', 'V13', 'V19', 'V25',...
        'V22', 'V23', 'mVUM301', 'V11', 'V08', 'V09', 'V05', 'V04', ...
        'V35', 'V43', 'V44', 'V07', 'V40'} ;
end
if ~exist('conn_clust_struct','var') || isempty(conn_clust_struct)
    % get path to current mfile so we can find clustering results by
    % relative location
    [mfilePath, ~, ~] = fileparts(mfilename('fullpath')) ;
    figDirectory = fileparts(mfilePath) ;
    parentDirectory = fileparts(figDirectory) ;
    clustPath = fullfile(parentDirectory, 'data', 'clustering') ;
    
    % switch depending on which neurons we're looking at
    if useINClusterFlag && strcmpi(neuronType,'IN')
        conn_clust_struct = importdata(fullfile(clustPath, ...
            'conn_clust_struct_IN.mat')) ;
    elseif useDNClusterFlag && strcmpi(neuronType,'DN')
        conn_clust_struct = importdata(fullfile(clustPath, ...
            'conn_clust_struct_DN.mat')) ;
    end
end

% -----------------------------------------------------------------
%% how we sort neurons depends on which neurons we're dealing with
switch neuronType
    case 'IN'
        % --------------------------------------------------------------
        % here we're dealing with interneurons. either use connectivity
        % clusters or sort alphabetically
        if useINClusterFlag
            % -------------------------------------------------
            % in this case, use IN clusters
            [~, ~, outperm] = get_cluster_outperm(conn_clust_struct) ;
            
            % read out labels from conn_clust_struct
            cluster_IN_names = {conn_clust_struct.neuron_label}' ;
            cluster_IN_names = cellfun(@(y) strrep(y,'_','-'), ...
                cluster_IN_names, 'UniformOutput',false) ;
            
            % first try to make the label_cell_1 and cluster_IN_names cells match
            % up
            [tf, cluster_sort_idx] = ismember(cluster_IN_names, labels) ;
            
            if sum(~tf) > 0
                fprintf('Error finding matching names! \n')
                keyboard
            end
            
            % now use cluster_sort_idx and outperm to make label_cell_1 match
            % dendrogram plot order
            sort_ind = cluster_sort_idx(outperm) ;
            labels_sort = labels(sort_ind) ;
        else
            % ---------------------------------------
            % otherwise just use sort
            [labels_sort, sort_ind] = sort(labels) ;
        end
        
    case {'MN', 'VUM'}
        % -----------------------------------------------
        % here we sort MNs or VUMs by pre-defined order
        
        % figure out whether we're dealing with MN or VUM
        if strcmp(neuronType, 'MN')
            neuronOrder = MN_order ;
        elseif strcmp(neuronType, 'VUM')
            neuronOrder = VUM_order ;
        else
            fprintf('Error: invalid neuron type: %s \n', neuronType)
            keyboard
        end
        
        % sort by list
        sort_ind = zeros(length(neuronOrder),1) ;
        for q = 1:length(neuronOrder)
            match_ind = find(cellfun(@(y) contains(y,neuronOrder{q},...
                'IgnoreCase',1), labels));
            if ~isempty(match_ind)
                sort_ind(q) = match_ind ;
            else
                fprintf('Error finding matching names! \n')
                keyboard
            end
        end
        labels_sort = labels(sort_ind) ;
        
        
    case 'DN'
        % --------------------------------------------------------------
        % here dealing with DESCENDING NEURONS. either use connectivity
        % clusters or sort alphabetically
        if useDNClusterFlag
            % -------------------------------------------------
            % in this case, use IN clusters
            [~, ~, outperm] = get_cluster_outperm(conn_clust_struct) ;
            
            % read out labels from conn_clust_struct
            cluster_DN_names = {conn_clust_struct.neuron_label}' ;
            cluster_DN_names = cellfun(@(y) strrep(y,'_','-'), ...
                cluster_DN_names, 'UniformOutput',false) ;
            
            % first try to make the label_cell_1 and cluster_IN_names cells match
            % up
            [tf, cluster_sort_idx] = ismember(cluster_DN_names, labels) ;
            
            if sum(~tf) > 0
                fprintf('Error finding matching names! \n')
                keyboard
            end
            
            % now use cluster_sort_idx and outperm to make label_cell_1 match
            % dendrogram plot order
            sort_ind = cluster_sort_idx(outperm) ;
            
            % need to add in indices for DNs that don't have i/o
            % annotations
            [~, non_clustered_ind] = setdiff(labels, cluster_DN_names) ;
            sort_ind = [sort_ind ; non_clustered_ind] ;
            
            % sort labels
            labels_sort = labels(sort_ind) ;
        else
            % ---------------------------------------
            % otherwise just use sort
            [labels_sort, sort_ind] = sort(labels) ;
        end
end

end