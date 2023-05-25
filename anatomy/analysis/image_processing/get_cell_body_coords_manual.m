% -------------------------------------------------------------------------
% function to get the position of a cell body from an image by clicking on
% coordinates
% -------------------------------------------------------------------------
function cell_body_coords = get_cell_body_coords_manual(bwMat) 
% create mip to make it easy to click on cell body
bwMIP = max(bwMat,[],3) ; 
h_temp = figure ; 
imshow(bwMIP) 


end