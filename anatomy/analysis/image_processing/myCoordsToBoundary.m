% -------------------------------------------------------------------------
% function to convert coordinates (indices of all "on" pixels in bw image)
% into a list of 3D coordinates and a boundary
%
% NB: we assume here that coords is just a list of points. if we load a
% mask coords file, e.g., we'll get a cell array (input and output), so
% we'd need to grab one of those as input here
% -------------------------------------------------------------------------
function [vox_coord_cell, vox_boundary_cell] = ...
    myCoordsToBoundary(coords, imSize, boundaryConvexity, debugFlag) 
% --------------------
%% inputs/params
if ~exist('imSize','var') || isempty(imSize)
   imSize = [1119, 573, 333] ; % this is the default image size 
end
if ~exist('boundaryConvexity','var') || isempty(boundaryConvexity)
   boundaryConvexity = 0.9 ; % 1.0 makes tightest boundary ; 0.0 makes convex hull
end
if ~exist('debugFlag','var') || isempty(debugFlag)
   debugFlag = false ; % check output?
end

cleanFlag = true ; % clean BW image?
smoothFactor = 0.1 ; % use smooth3 to apply some gaussian smoothing to perim
% -----------------------------------------------
%% generate fake image and get 3D coords (sub)
% NB: normally we could just use ind2sub, but we want to just get the
% perimeter here. 

% generate bw image
bwAll = false(imSize) ; 
bwAll(coords) = true ; 

% get number of connected components in image
CC = bwconncomp(bwAll) ; 
N_obj = CC.NumObjects ; 

% -------------------------------------------------------
%% get boundary info from images
% just use the function we already wrote...

% initialize output
vox_coord_cell = cell(N_obj,1) ;
vox_boundary_cell = cell(N_obj,1) ;

% perform loop over connected components
for k = 1:N_obj
    % generate bw image of just this cc
    bw = false(imSize) ;
    bw(CC.PixelIdxList{k}) = true ;
    
    % get coords and boundary
    [vox_coord_cell{k}, vox_boundary_cell{k}] = getVoxelBoundaryPoints(bw, ...
        smoothFactor, boundaryConvexity, cleanFlag, debugFlag) ;
end

end