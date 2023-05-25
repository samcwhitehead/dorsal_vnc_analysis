% -------------------------------------------------------------------------
% function to draw a neuron to existing axis. should hopefully be generally
% useful
% -------------------------------------------------------------------------
function ax = drawNeuronToAxis(ax, neuronData, neuronColor, neuronAlpha, ...
    imSize, parent)
% ------------------------------------
%% inputs
if ~exist('ax', 'var') || isempty(ax)
    ax = gca ;
end
if ~exist('neuronColor', 'var') || isempty(neuronColor)
    neuronColor = 0.5*[1.0, 1.0, 1.0] ; % default to black
end
if ~exist('neuronAlpha', 'var') || isempty(neuronAlpha)
    neuronAlpha = 0.5; % default to half transparent
end
% general image size -- will be used for converting coords to image
if ~exist('imSize', 'var') || isempty(imSize)
    imSize = [1119, 573, 333]; % this is the default image size
end
if ~exist('parent', 'var') || isempty(parent)
    parent = [] ; 
end

% current axis properties
holdFlag = ishold(ax) ;
hold(ax,'on')

% dorsal view settings
az_dorsal = -180 ; 
el_dorsal = 90 ;

% level to use for isosurface
iso_level = 0.5 ; % for a binary image, 0.5 should be right on edge (i.e. between "1" (true) and "0" (false))
% ----------------------------------------------------
%% determine what type of data we're dealing with
% first check if this is a skeleton structure
if isstruct(neuronData)
    if isfield(neuronData, 'skel')
        fprintf(['Input is a skeleton struct -- should just use'  ...
            '"drawSkelGraph.m ... \n'])
        ax = drawSkelGraph(ax, neuronData, neuronColor, [], neuronAlpha) ;
        return
    else
        fprintf('Error: unrecognized input type \n')
        keyboard
    end
else
    % if we get here, our data is some form of array. so we should now just
    % check whether it's a list of 3D coordinates or a binary image
    dataSize = size(neuronData) ;
    
    % if the data has 3 array dimensions, assume it's an image. if it has 2
    % assume it's a list of 3D coordinates (in subscript form). if it has 1
    % assume it's a list of 3D coordinates (in index form).
    if length(dataSize) == 3
        % in this case, we have a bw image (probably). just check that it's
        % logical and resize if need be
        if ~islogical(neuronData)
            fprintf('Error: was expecting binary image for 3D input \n')
            keyboard
        end
        if ~all(dataSize == imSize)
            neuronData = imresize3(neuronData, imSize) ;
        end
        
        % just explicitly name it BW to match up with the rest of the code
        neuronBW = neuronData ;
        
    elseif (length(dataSize) == 2) && all(dataSize > 1)
        % in this case, we assume we have a list of xyz coordinates.
        % convert these to index and generate a bw image
        neuronBW = false(imSize) ;
        bw_idx = sub2ind(imSize, neuronData(:,2), neuronData(:,1), ...
            neuronData(:,3)) ;
        neuronBW(bw_idx) = true ;
        
    elseif length(dataSize) == 2 && ~all(dataSize > 1)
        % in this case, we assume we have a set of coordinate indices. use
        % this to generate a bw image
        neuronBW = false(imSize) ;
        neuronBW(neuronData) = true ;
    else
        fprintf('Error: weird input \n')
        keyboard
        
    end
end

% --------------------------------------------------------------
%% draw neuron based on BW image
% if we've gotten  here, we should have a binary image corresponding to our
% neuron. calculate the isosurface of this bw neuron image and draw
neuronIsoSurf = isosurface(neuronBW, iso_level) ;

% plot neuron isosurface
neuronPatch = patch(ax, neuronIsoSurf) ;

% set neuron patch propoerties
neuronPatch.FaceColor = neuronColor ;
neuronPatch.FaceAlpha = neuronAlpha ;
neuronPatch.LineStyle = 'none' ;
%     neuronPatch.EdgeColor = colorCurr ;
%     neuronPatch.FaceAlpha = neuronAlpha ;
%     neuronPatch.EdgeAlpha = neuronAlpha ;

% ---------------------------------------------------------------
%% set some axis properties (but not too many)
axis(ax,'equal') ; % same units across x,y,z

% set view (dorsal by default)
view(az_dorsal,el_dorsal) ;

% try to ensure that things don't change when we alter view
set(ax,'DataAspectRatio',[1 1 1])
set(ax,'PlotBoxAspectRatio',[1 1 1])
camva('manual')

% remove axis lines/ticks/labels
axis off

% return axis to previous hold state
if ~holdFlag
   hold(ax, 'off') 
end

if ~isempty(parent)
    grp = hgtransform('Parent',parent) ;
    set(neuronPatch,'Parent',grp) ;  
end
% assume we'll want to rotate
% rotate3d on
end