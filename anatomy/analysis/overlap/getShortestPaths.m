% -------------------------------------------------------------------------
% quick and dirty script to find shortest paths between neurons of
% different sets
%
% NB: need to generate combined overlap matrices for this
% (combineOverlapMats.m)
% -------------------------------------------------------------------------
%% path info and params
[mfilePath, ~, ~] = fileparts(mfilename('fullpath')) ;
figDirectory = fileparts(mfilePath) ;
parentDirectory = fileparts(figDirectory) ;
rootPath = fullfile(parentDirectory, 'data') ;
combOverlapPath = fullfile(rootPath, 'overlap_calculations', 'masked', ...
    'combined') ;
graphFn = 'overlap_digraph_all.mat' ; 

% source and destination neuron types
sourceNeuronType = 'DN' ; 
destNeuronType = 'MN' ; 

% binarize graph so that we remove edges? if not, need to invert weights,
% since the default is to have higher weight count as more distance
pathMethod = 'unweighted' ; % 'positive' or 'unweighted' 

% save output?
saveFlag = true ;
savePath = combOverlapPath ; 

% -----------------------------------------------------------------------
%% load overlap graph and neuron type struct
overlap_graph = importdata(fullfile(combOverlapPath, graphFn)) ; 
neuronTypeStruct = importdata(fullfile(rootPath, 'neuronTypeStruct.mat')) ; 

% get node names in struct field form
node_names = overlap_graph.Nodes.Name ; 
node_names_sf = cellfun(@(y) strrep(y,'-','_'), node_names, ...
    'UniformOutput',false) ; 

% find indices for source and destination neuron types
source_idx = cellfun(@(y) strcmpi(neuronTypeStruct.(y), ...
    sourceNeuronType), node_names_sf) ; 
dest_idx = cellfun(@(y) strcmpi(neuronTypeStruct.(y), ...
    destNeuronType), node_names_sf) ; 

% ------------------------------------------------------------------------
%% binarize graph edges?
% read out current edge weights (should correspond to overlap voxel count)
edge_weights = overlap_graph.Edges.Weight ; 

% if binarizeFlag
%     % if binarizing, just convert to 0 or 1
%     overlap_graph.Edges.Weight = double(edge_weights > 0) ; 
if strcmp(pathMethod, 'positive')
   % otherwise, take inverse of weights so that more overlap => less distance
   non_zero_idx = (edge_weights > 0) ; 
   edge_weights(non_zero_idx) = (edge_weights(non_zero_idx)).^(-1) ; 
   overlap_graph.Edges.Weight = edge_weights ; 
end

% -----------------------------------------------------------------------
%% loop over source and destination options
% ...find shortest paths, put them in a table (each entry will be a CELL)
% initialize storage
Ns = sum(source_idx) ; 
Nd = sum(dest_idx) ; 

shortest_path_table = table('Size',[Ns, Nd], 'VariableTypes', ...
    repmat(["string"],1,Nd)) ;
shortest_path_table.Properties.RowNames = node_names(source_idx) ; 
shortest_path_table.Properties.VariableNames = node_names(dest_idx) ; 

% get source and dest idx in subscript form
source_ind = find(source_idx) ; 
dest_ind = find(dest_idx) ; 

% begin loop over source neurons
for ii = 1:Ns
    % get info of source neuron 
    sourceNeuronID = node_names(source_ind(ii)) ; 
    sourceNeuronNum = findnode(overlap_graph, sourceNeuronID) ; 
    
    % loop over destination neurons
    for jj = 1:Nd
        % info for current destination neuron
        destNeuronID = node_names(dest_ind(jj)) ; 
        destNeuronNum = findnode(overlap_graph, destNeuronID) ; 
        
        % find shortest path between current source/destination pair
        [P, d] = shortestpath(overlap_graph, sourceNeuronNum, ...
            destNeuronNum, 'Method', pathMethod) ; 
        
        % add info to  table
        shortest_path_table{ii,jj} = string(strjoin({node_names{P}},', ')) ;  % strjoin({node_names{P}},', ') ; 
        
    end
end

% -----------------------------------------------------------
%% save results?
if saveFlag
    saveFn = sprintf('shortest_path_%s_%s_%s.mat', sourceNeuronType, ...
        destNeuronType, pathMethod) ;
    save(fullfile(savePath, saveFn),'shortest_path_table')
end

