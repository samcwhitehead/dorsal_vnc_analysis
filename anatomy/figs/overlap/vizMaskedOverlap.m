% -------------------------------------------------------------------------
% script to make figure showing neuron overlap with input/output masks
% visualized as well
% -------------------------------------------------------------------------
%% plot params
figPosition = [5.4583, 1.8333, 6.2917, 7.8125] ; %[7.0000, 5.7292, 1.06, 3.14] ;
figUnits = 'inches' ;

clusterColors =  brewermap(15, 'Dark2') ;
inputColor = clusterColors(1,:) ; %  %
outputColor = [0, 0, 0] ;
overlapColor = [0.8, 0, 0] ;
maskColor = 0.7*[1, 1, 1] ;

neuronAlpha = 0.85 ;
overlapAlpha = 1.0 ; 
% notInMaskAlpha = 0.20 ;
maskAlpha = 0.15 ;
vncAlpha = 0.00 ;

% general params
drawOverlapFlag = true ; % color overlapping regions separately?
drawVNCFlag = true ;  % draw VNC mask regions?

INPUT = 1 ;
OUTPUT = 2 ;

imSize = [1119, 573, 333] ; % this is the default image size

% save info
saveFlag = false ;
[mfilePath, ~, ~] = fileparts(mfilename('fullpath')) ; 
savePath = fullfile(mfilePath, 'output') ;

% ----------------------------------
%% load neuron data
% cell arrays containing 1) neuron name and 2) neuron type
input_neuron =  {'WBL008', 'IN'} ; % {'WBL014', 'IN'} ; % {'XBL01', 'IN'} ; %
output_neuron = {'DNp01', 'DN'} ; % {'DNg02', 'DN'} ; % 

dataType = 'coords_sym' ; % use coords or bw image? symmetrized or no?
maskType = 'mask_coords_sym' ; % same question for i/o masks

% load data for neurons
fprintf('Loading neuron and mask data ... \n')

input_coords = myLoadNeuronData(input_neuron{1}, 'coords_non_sym', input_neuron{2});
output_coords = myLoadNeuronData(output_neuron{1},dataType,output_neuron{2});

% load data for masks
input_masks = myLoadNeuronData(input_neuron{1}, maskType, input_neuron{2});
output_masks = myLoadNeuronData(output_neuron{1},maskType,output_neuron{2});

% ----------------------------------------------------------
%% process data to get the things we want to plot
fprintf('Processing neuron and mask data ... \n')

% combine input mask of input neuorn and output mask of output neuron:
mask_comb = intersect(input_masks{INPUT}, output_masks{OUTPUT}) ;

% get portion of both neurons that are within mask
if drawOverlapFlag
    input_in_mask = intersect(input_coords, mask_comb) ;
    output_in_mask = intersect(output_coords, mask_comb) ;
    
    % % ... as well as sections outside mask
    % input_not_in_mask = setdiff(input_coords, input_in_mask) ;
    % output_not_in_mask = setdiff(output_coords, output_in_mask) ;
    
    % use the above to get coordinates where images actually overlap
    % (within mask)
    overlap_coords = intersect(input_in_mask, output_in_mask) ;
end

% finally, get boundary points for mask neuropil regions
[mask_coord_cell, mask_boundary_cell] = myCoordsToBoundary(mask_comb) ;

% ----------------------------------------------------------------------
%% start to plot stuff
fprintf('Generating figure ... \n')

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
        vncAlpha, parent)
end

% % ---------------------------------
% % draw neurons outside mask
% ax = drawNeuronToAxis(ax, input_not_in_mask, inputColor, notInMaskAlpha) ;
% ax = drawNeuronToAxis(ax, output_not_in_mask, outputColor, notInMaskAlpha) ;

% ---------------------------------
% draw neurons
ax = drawNeuronToAxis(ax, input_coords, inputColor, neuronAlpha, ...
    imSize, parent) ;
ax = drawNeuronToAxis(ax, output_coords, outputColor, neuronAlpha, ...
    imSize, parent) ;

% ---------------------------------
% draw overlap sections?
if drawOverlapFlag
    ax = drawNeuronToAxis(ax, overlap_coords, overlapColor, overlapAlpha, ...
        imSize, parent) ;
end

% --------------------------------
% draw neuropil mask region
N_mask_obj = length(mask_boundary_cell) ;
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

% center objects then set axis limits
T = makehgtform('translate', -0.5*[imSize(2), imSize(1), imSize(3)]) ;
set(parent, 'Matrix', T)
set(ax, 'xlim', 0.5*imSize(2).*[-1,1], 'ylim', 0.5*imSize(1).*[-1,1], ...
    'zlim', 0.5*imSize(3).*[-1,1])

% add lighting?
camlight;
lighting phong
% material metal

fprintf('Completed drawing neurons/mask \n')

% --------------------------------------------------
% save figure
if saveFlag
    fprintf('Saving figure ... \n')
    % base savename
    saveFn = ['ex_overlap_' input_neuron{1} '_' output_neuron{1}] ;
    if drawVNCFlag
        if vncAlpha > 0.0
            saveFn = [saveFn '_withVNC'] ;
        else
            saveFn = [saveFn '_withTransparentVNC'] ;
        end
        
    end
    if drawOverlapFlag
       saveFn = [saveFn '_withOverlap'] ;  
    end
    % save dorsal view
    export_fig(h_main, fullfile(savePath, [saveFn '_dorsal']), ...
        '-dpng','-r300')
    savefig(h_main, fullfile(savePath, [saveFn '.fig'])) ;
    
    % change view to sagittal and save
    Ry = makehgtform('yrotate', pi/2) ;
    set(parent, 'Matrix', Ry*T)
    export_fig(h_main, fullfile(savePath, [saveFn '_sagittal']), ...
        '-dpng','-r300')
    
    % change view to transverse and save
    Rx = makehgtform('xrotate', pi/2) ;
    %Rz = makehgtform('zrotate', pi) ;
    set(parent, 'Matrix', Rx*T)
    export_fig(h_main, fullfile(savePath, [saveFn '_transverse']), ...
        '-dpng','-r300')
    
    % change view to ventral
    Ry = makehgtform('yrotate', pi) ;
    set(parent, 'Matrix', Ry*T)
    export_fig(h_main, fullfile(savePath, [saveFn '_ventral']), ...
        '-dpng','-r300')
end
