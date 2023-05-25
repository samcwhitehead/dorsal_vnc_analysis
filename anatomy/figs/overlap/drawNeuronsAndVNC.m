% -------------------------------------------------------------------------
% function to generate image of neuron(s) with VNC as background
% -------------------------------------------------------------------------
function ax = drawNeuronsAndVNC(ax, fig, vncBW, neuronBWCell, ...
    neuronColors, vncColor, vncAlpha)
% -------------------------------------
%% inputs and params
if ~exist('ax','var') || isempty(ax)
   ax = gca ;  
end
if ~exist('fig','var') || isempty(fig)
   fig = gcf ;  
end
if ~exist('neuronColors','var') || isempty(neuronColors)
   neuronColors = zeros(length(neuronBWCell), 3) ;  
end
if ~exist('vncColor','var') || isempty(vncColor)
    vncColor = 0.7*[1, 1, 1] ; 
end
if ~exist('vncAlpha','var') || isempty(vncAlpha)
    vncAlpha = 0.05 ;
end

% general plot params
az_dorsal = -180 ; 
el_dorsal = 90 ;

neuronAlpha = 1.0 ; 

shrinkFactor = 7 ; % how much to reduce volume by 
% -------------------------------------------------------------------------
%% process VNC (get perimeter, shrink, get coordinates + boundary)
if ~isempty(vncBW)
    vncPerim = bwperim(vncBW) ; 
    vncPerim_smooth = imbinarize(smooth3(vncPerim, 'gaussian',9)) ;
    nv = reducevolume(vncPerim_smooth, shrinkFactor) ; 
    
    [qy, qx, qz] = ind2sub(size(nv), find(nv(:))) ;
    % image <-> cartesian
    %qy = size(nv,1) - qy + 1 ;
    %qz = size(nv,3) - qz + 1 ;
    vnc_bound = boundary(qx, qy, qz) ;
else
    qx = nan ; 
    qy = nan ; 
    qz = nan ; 
    vnc_bound = nan ; 
end

% -------------------------------------------------------------------------
%% process neurons (get isosurfaces
N_neurons = length(neuronBWCell) ; 
neuronIsoCell = cell(N_neurons, 1) ; 
for i = 1:N_neurons
    % calculate neuron isosurface (0.5 is the boundary between 0 and 1 for
    % bw images
    neuronBW = neuronBWCell{i} ;
    fv = isosurface(neuronBW,0.5) ;
    neuronIsoCell{i} = fv ; 
end

% -------------------------------------------------------------------------
%% draw objects
% get current fig/axis properties
set(fig,'CurrentAxes',ax)
holdFlag = ishold(ax) ; 
hold(ax,'on') 

% ----------------------------
% draw vnc (triangulation)
trisurf(vnc_bound, qx*shrinkFactor, qy*shrinkFactor, qz*shrinkFactor,...
    'FaceColor',vncColor, 'EdgeColor','none','FaceAlpha',vncAlpha) ; 
    
% ------------------------------
% draw neurons
for j = 1:N_neurons
    % current plot color
    colorCurr = neuronColors(j,:) ; 
    
    % plot neuron isosurface
    p_neuron = patch(neuronIsoCell{j}) ;
    
    % set neuron patch propoerties
    p_neuron.FaceColor = colorCurr ;
    p_neuron.EdgeColor = colorCurr ;
    p_neuron.FaceAlpha = neuronAlpha ;
    p_neuron.EdgeAlpha = neuronAlpha ;
      
end

% axis properties
axis(ax, 'equal')
camlight
lighting gouraud
view(az_dorsal,el_dorsal)
axis off

if holdFlag
    hold(ax ,'on')
else
    hold(ax, 'off')
end
% message = sprintf([line_name '\n' driver]) ;
% text(max(shrinkFactor*qx)+20, min(shrinkFactor*qy) - 20,...
%     max(shrinkFactor*qz), message,'fontSize',16)

end