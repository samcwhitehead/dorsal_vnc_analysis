% -------------------------------------------------------------------------
% function to visualize a path, i.e. a set of neurons connecting a selected
% source and destination neuron
%{
% EXAMPLE USAGE:
sourceNeuron = 'DNg02' ; 
destNeuron = 'DLMn5' ; 

sourceNeuronType = 'DN' ; 
destNeuronType = 'MN' ; 

rootPath = 'D:\Fly Imaging\Erica Interneuron Stacks\' ; % path to data
combOverlapPath = fullfile(rootPath, 'overlap_calculations', 'masked', ...
    'combined') ;
pathMethod =  'positive'  ; % 'positive' | 'unweighted'

% load data 
dataFn = sprintf('shortest_path_%s_%s_%s.mat', sourceNeuronType, ...
        destNeuronType, pathMethod) ;
path_table = importdata(fullfile(combOverlapPath, dataFn)) ; 

[h_main, ax, parent] = visualizeNeuronPath(sourceNeuron,destNeuron, path_table) ;
%}
% -------------------------------------------------------------------------
function [h_main, ax, parent] = visualizeNeuronPath(sourceNeuron, ...
    destNeuron, path_table, cmap)
% ------------------------------
%% 
% inputs and params
if ~exist('path_table','var') || isempty(path_table)
    % table containing shortest distances between sets of neurons
   path_table_path = ['D:\Fly Imaging\Erica Interneuron Stacks\' ...
       'overlap_calculations\masked\combined\'] ; 
   path_table_fn = 'shortest_path_DN_MN_unweighted.mat' ; 
   path_table = importdata(fullfile(path_table_path, path_table_fn)) ; 
end
if ~exist('cmap','var') || isempty(cmap)
   cmap = lines(10) ;  
end

% delimiter to read out paths from table
entry_delim = ', ' ; 

% neuron data type
dataType = 'coords_non_sym' ; 

% image size
imSize = [1119, 573, 333] ; % this is the default image size

% some plot params
neuronAlpha = 0.25 ;

% ---------------------------------------------
%%
% find indices in path_table for source and destination neurons
fprintf('Getting neuron path ... \n')
sourceNames = path_table.Properties.RowNames ; 
destNames = path_table.Properties.VariableNames ; 

source_idx = cellfun(@(y) strcmpi(sourceNeuron, y), sourceNames) ; 
dest_idx = cellfun(@(y) strcmpi(destNeuron, y), destNames) ; 

% get path from selected source and destination
path_curr = path_table{source_idx, dest_idx} ; 
path_split = strsplit(path_curr, entry_delim) ; 
N_path = length(path_split) ; % number of neurons in path (including start and end)

% ----------------------------------------------
%%
% Load data for neurons along path
fprintf('Loading neuron data ... \n')
neuronDataCell = cell(N_path,1) ; 
for k = 1:N_path
    neuronDataCell{k} = myLoadNeuronData(path_split{k}, dataType) ; 
    if isempty(neuronDataCell{k})
        dataTypeTemp = strrep(dataType,'non_sym','sym') ; 
        neuronDataCell{k} = myLoadNeuronData(path_split{k}, dataTypeTemp) ;
    end
end

% ------------------------------------------------
%%
% Generate figure
fprintf('Making figure ... \n')
h_main = figure ; 

% initialize axis
ax = gca ;
hold(ax,'on')

% create a parent hgtform to hopefully use to move stuff around later
parent = hgtransform('Parent',ax);

% ---------------------------------
% draw VNC as background
drawVNCBackgroundToAxis(ax, h_main, [], [], [], parent) ; 

% draw neurons from path 
for m = 1:N_path
    ax = drawNeuronToAxis(ax, neuronDataCell{m}, cmap(m,:), neuronAlpha, ...
        imSize, parent) ;
end

% ------------------------------------
%%
% set some axis properties
set(ax,'DataAspectRatio',[1 1 1])
set(ax,'PlotBoxAspectRatio',[1 1 1])
camva('manual')

% center objects then set axis limits
T = makehgtform('translate', -0.5*[imSize(2), imSize(1), imSize(3)]) ;
set(parent, 'Matrix', T)
set(ax, 'xlim', 0.5*imSize(2).*[-1,1], 'ylim', 0.5*imSize(1).*[-1,1], ...
    'zlim', 0.5*imSize(3).*[-1,1])

% change to sagittal view
Ry = makehgtform('yrotate', pi/2) ;
set(parent, 'Matrix', Ry*T)
    
% add lighting?
camlight;
lighting phong
% material metal

fprintf('Completed drawing neurons \n')

end