% -------------------------------------------------------------------------
% function to get bw image from structure containing a boundary and coords
% (output from "getVoxelBoundaryPoints.m")
% -------------------------------------------------------------------------
function bw_struct = getBWfromCoords(bw_struct, imDim) 
% inputs
if ~exist('imDim','var') || isempty(imDim)
    imDim = [] ; % will fill in later
end

% read out coords from struct
coords = bw_struct.vox_coord_cell{:} ;

% get image dimensions from coords if we don't have them as input
if isempty(imDim)
   max_x = max(coords(:,1)) ; 
   max_y = max(coords(:,2)) ; 
   max_z = max(coords(:,3)) ; 
   imDim = [max_y, max_x, max_z] ; 
end

% initialize binary image
BW = false(imDim) ; 

% get indices for coordinate locations -- note the 1<->2 dim switch
idx = sub2ind(imDim, coords(:,2), coords(:,1), coords(:,3)) ; 

% fill image at coords points
BW(idx) = true ; 

% % note that these are just the boundary coordinates, so try to fill image as
% % well
% BW = bwmorph3(BW,'fill') ; 

% add image to struct
bw_struct.BW = BW ; 
end