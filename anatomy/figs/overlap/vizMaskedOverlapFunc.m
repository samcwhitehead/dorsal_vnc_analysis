% -------------------------------------------------------------------------
% FUNCTION to make figure showing neuron overlap with input/output masks
% visualized as well
% -------------------------------------------------------------------------
function h_main = vizMaskedOverlapFunc(input_neuron, output_neuron, ...
    inputColor, outputColor, vncAlpha, drawVNCFlag, drawMaskFlag, ...
    drawOverlapFlag, labelFlag, saveFlag, savePath, figPosition)
% -----------------------------------------
%% inputs and plot params
if ~exist('input_neuron','var') || isempty(input_neuron)
    input_neuron =  {'WBL13', 'IN', 'coords_non_sym'} ; %{'WBL14', 'IN'} ; % {'XBL01', 'IN'} ; %
end
if ~exist('output_neuron','var') || isempty(output_neuron)
    output_neuron = {'DNp01', 'DN', 'coords_sym'} ;
end
if ~exist('inputColor','var') || isempty(inputColor)
    clusterColors =  brewermap(15, 'Dark2') ;
    inputColor = clusterColors(1,:) ; %  %
end
if ~exist('outputColor','var') || isempty(outputColor)
    outputColor = [0, 0, 0] ;
end
if ~exist('vncAlpha','var') || isempty(vncAlpha)
    vncAlpha = 0.05 ;
end
if ~exist('drawVNCFlag','var') || isempty(drawVNCFlag)
    drawVNCFlag = true ;
end
if ~exist('drawMaskFlag','var') || isempty(drawMaskFlag)
    drawMaskFlag = true ;
end
if ~exist('drawOverlapFlag','var') || isempty(drawOverlapFlag)
    drawOverlapFlag = false ; % color overlapping regions separately?
end
% label image with neuron names?
if ~exist('labelFlag','var') || isempty(labelFlag)
    labelFlag = true ;
end
% save info
if ~exist('saveFlag','var') || isempty(saveFlag)
    saveFlag = false ;
end
if ~exist('savePath','var') || isempty(savePath)
    savePath = ['D:\Dropbox\Paper Manuscripts\Janelia vnc paper\' ...
        'figure drafts\DN_overlap_fig\overlap_illustration\'];
end
% fig outer position
if ~exist('figPosition','var') || isempty(figPosition)
    figPosition = [5.4583, 1.8333, 6.2917, 7.8125] ; 
end

% --------------------------------------------------
% some plot params
figUnits = 'inches' ;

overlapColor = [0.8, 0, 0] ;
maskColor = 0.7*[1, 1, 1] ;

neuronAlpha = 0.85 ;
% notInMaskAlpha = 0.20 ;
maskAlpha = 0.15 ;

% which data type to use for i/o masks
maskType = 'mask_coords_sym' ;

% general params
INPUT = 1 ;
OUTPUT = 2 ;

imSize = [1119, 573, 333] ; % this is the default image size

labelFontSize = 6 ; % assuming we can get to a point where we don't need to scale drawings in Inkscape...
% ----------------------------------
%% load neuron data
% load data for neurons
fprintf('Loading neuron and mask data ... \n')

input_coords = myLoadNeuronData(input_neuron{1}, input_neuron{3}, ...
    input_neuron{2});
output_coords = myLoadNeuronData(output_neuron{1}, output_neuron{3},...
    output_neuron{2});

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
        vncAlpha, parent) ;
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
    ax = drawNeuronToAxis(ax, overlap_coords, overlapColor, neuronAlpha, ...
        imSize, parent) ;
end

% --------------------------------
% draw neuropil mask region
if drawMaskFlag
    N_mask_obj = length(mask_boundary_cell) ;
    grp_array = gobjects(N_mask_obj,1) ;
    
    for k = 1:N_mask_obj
        grp_array(k) = draw3Dboundary(ax, parent, mask_coord_cell{k}, ...
            mask_boundary_cell{k}, maskColor, maskAlpha) ;
    end
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

% add neuron names as label?
if labelFlag
    % coordinates for label (in units of data)
    text_x = -0.4*imSize(2) ; %-0.35*imSize(2) ;
    text_y = -0.425*imSize(1) ;
    text_z = 0.425*imSize(3) ;
    
    % string to display
    txt_str = sprintf(['{{\\color[rgb]{%f, %f, %f} %s}}\n' ...
        '{{\\color[rgb]{%f, %f, %f} %s}}'], ...
        inputColor(1), inputColor(2), inputColor(3), input_neuron{1},...
        outputColor(1), outputColor(2), outputColor(3), output_neuron{1}) ;
    
    % add text to axis
    txt = text(ax, text_x, text_y, text_z, txt_str, ...
        'FontSize', labelFontSize, ...
        'FontName','arial',...
        'HorizontalAlignment','right') ;
end
fprintf('Completed drawing neurons/mask \n')

% --------------------------------------------------
% save figure
if saveFlag
    fprintf('Saving figure ... \n')
    % -----------------------------------------------
    % base savename
    saveFn = [input_neuron{1} '_' output_neuron{1}] ;
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
    
    % ----------------------------------------------------
    % save dorsal view
    exportgraphics(h_main, fullfile(savePath, [saveFn '_dorsal.png']), ...
        'Resolution',500)
    savefig(h_main, fullfile(savePath, [saveFn '.fig'])) ;
    
    % ------------------------------------------------------
    % change view to sagittal and save
    Ry = makehgtform('yrotate', pi/2) ;
    set(parent, 'Matrix', Ry*T)
    
    % adjust label, if using
    if labelFlag
        % delete old label
        delete(txt)
        
        % add text to axis, switching coordinates around
        txt = text(ax, 0.5*text_x, text_y, text_z, txt_str, ...
            'FontSize', labelFontSize, ...
            'FontName','arial',...
            'HorizontalAlignment','right') ;
    end
    exportgraphics(h_main, fullfile(savePath, [saveFn '_sagittal.png']), ...
        'Resolution', 500)
    
    % --------------------------------------------------------
    % change view to transverse and save
    Rx = makehgtform('xrotate', pi/2) ;
    %Rz = makehgtform('zrotate', pi) ;
    set(parent, 'Matrix', Rx*T)
    
   set(ax, 'xlim', 0.5*imSize(2).*[-1,1], 'zlim', 0.5*imSize(1).*[-1,1], ...
    'ylim', 0.5*imSize(3).*[-1,1])

    % adjust label, if using
    if labelFlag
        % delete old label
        delete(txt)
        
        % add text to axis, switching coordinates around
        txt = text(ax, 0.9*text_x, -1.5*text_z, 0, txt_str, ...
            'FontSize', labelFontSize, ...
            'FontName','arial',...
            'HorizontalAlignment','right') ;
    end
    
    exportgraphics(h_main, fullfile(savePath, [saveFn '_transverse.png']), ...
        'Resolution', 500)
    
    % --------------------------------------------------------
    % change view to ventral and save
    Ry = makehgtform('yrotate', pi) ;
    set(parent, 'Matrix', Ry*T)
    
    set(ax, 'xlim', 0.5*imSize(2).*[-1,1], 'ylim', 0.5*imSize(1).*[-1,1], ...
    'zlim', 0.5*imSize(3).*[-1,1])

    % adjust label, if using
    if labelFlag
        % delete old label
        delete(txt)
        
        % add text to axis, switching coordinates around
        txt = text(ax, text_x, text_y, text_z, txt_str, ...
            'FontSize', labelFontSize, ...
            'FontName','arial',...
            'HorizontalAlignment','right') ;
    end
    
    exportgraphics(h_main, fullfile(savePath, [saveFn '_ventral.png']), ...
        'Resolution',500)
end
