% -------------------------------------------------------------------------
% quick and dirty function to get the overlap between a neuron and a given
% neuropil mask (i.e. how many voxels does the neuron have in a given
% neuropil)
%
% NB: this is just raw voxel count -- we don't care if it's cell body,
% dendrite, axon, etc.
%
%
%{
neuron_name = 'HBI013' ; 
neuropil_name = 'HTct' ; 

[overlap_val, overlap_val_norm] = neuronNeuropilOverlap(neuron_name, neuropil_name) ; 

%}
% -------------------------------------------------------------------------
function [overlap_val, overlap_val_norm] = ...
    neuronNeuropilOverlap(neuron_name, neuropil_name, rootPath, ...
        maskPath, plotFlag)
% ---------------------------------
%% inputs and params
if ~exist('rootPath','var') || isempty(rootPath)
    rootPath = 'D:\Fly Imaging\Erica Interneuron Stacks\' ;
end
if ~exist('maskPath','var') || isempty(maskPath)
    maskPath = fullfile(rootPath, ...
        'VNC neuropil','segmented_images','registered') ; 
end
if ~exist('plotFlag','var') || isempty(plotFlag)
    % visualize results?
    plotFlag = false ;
end

% which type of data to grab
neuronDataType = 'coords_sym' ; 
neuropilDataType = 'coords' ; 

% ------------------------------------------------------------
%% load neuron and mask data 
% load neuron
neuronData = myLoadNeuronData(neuron_name, neuronDataType) ; 

% load mask
if strcmp(neuropilDataType,'coords')
    maskFn = fullfile(maskPath, sprintf('%s_coords.mat', neuropil_name)) ;
else
    maskFn = fullfile(maskPath, sprintf('%s.mat', neuropil_name)) ;
end
maskData = importdata(maskFn) ; 

% -------------------------------------------------------
%% calculate overlap
overlap_val = 0 ; 
for k = 1:length(maskData)
    overlap_val = overlap_val + numel(intersect(neuronData, maskData{k})) ; 
end

% also get overlap normalized by total number of neuron voxels
overlap_val_norm = overlap_val/numel(neuronData) ; 


% ------------------------------------------------------------
%% visualize neuron and neuropil?
if plotFlag
    % should write helper function below to keep code readable
    % fprintf('Under construction! \n')
    [h_main, ax] = vizNeuronAndNeuropil(neuronData, maskData) ;
    title(ax, sprintf('%s vs %s \n frac. overlap = %f', neuron_name, ...
        neuropil_name, overlap_val_norm))
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HELPER FUNCTION(S)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------------------------------------------------------
%% function to plot both neuron and neuropil mask
function [h_main, ax] = vizNeuronAndNeuropil(neuronData, maskData, imSize)
% -----------------------
% inputs and params
if ~exist('imSize','var') || isempty(imSize)
   imSize = [1119, 573, 333] ; % this is the default image size 
end

% some plot params
figPosition = [5.4583, 1.8333, 6.2917, 7.8125] ; %[7.0000, 5.7292, 1.06, 3.14] ;
figUnits = 'inches' ;

drawVNCFlag = true ; 

neuronColor = 0.0*[1, 1, 1] ;
maskColor = 0.7*[1, 0, 0] ;
neuronAlpha = 0.85 ;
maskAlpha = 0.15 ;
vncAlpha = 0.05 ; 

% ---------------------------------------------
% convert mask data into boundary for plotting
[mask_coord_cell, mask_boundary_cell] = myCoordsToBoundary(maskData{:}) ;

% ------------------------------------
% make figure 
% initialize figure window
h_main = figure('PaperPositionMode','auto','MenuBar','none',...
    'ToolBar','none','DockControls','off','Units',figUnits,...
    'OuterPosition',figPosition) ;

% initialize axis
ax = gca ;
hold(ax,'on')

% create a parent hgtform to hopefully use to move stuff around later
parent = hgtransform('Parent',ax);

% ---------------------------------
% draw VNC as background
if drawVNCFlag
    drawVNCBackgroundToAxis(ax, h_main, [], [], ...
        vncAlpha, parent) ;
end

% ---------------------------------
% draw neurons
ax = drawNeuronToAxis(ax, neuronData, neuronColor, neuronAlpha, ...
    imSize, parent) ;

% --------------------------------
% draw neuropil masks
N_mask_obj = length(mask_coord_cell) ;
grp_array = gobjects(N_mask_obj,1) ;

for k = 1:N_mask_obj
    grp_array(k) = draw3Dboundary(ax, parent, mask_coord_cell{k}, ...
        mask_boundary_cell{k}, maskColor, maskAlpha) ;
end

% ------------------------------------
% some axis properties
set(ax,'DataAspectRatio',[1 1 1])
set(ax,'PlotBoxAspectRatio',[1 1 1])
camva('manual')

% ideally make bg transparent?
set(ax,'Color','none') ;

% center objects then set axis limits
T = makehgtform('translate', -0.5*[imSize(2), imSize(1), imSize(3)]) ;
set(parent, 'Matrix', T)
set(ax, 'xlim', 0.5*imSize(2).*[-1,1], 'ylim', 0.5*imSize(1).*[-1,1], ...
    'zlim', 0.5*imSize(3).*[-1,1])

% add lighting?
camlight;
lighting phong
% material metal

% --------------------------------------------------------
% change view to sagittal 
Ry = makehgtform('yrotate', pi/2) ;
set(parent, 'Matrix', Ry*T)

end