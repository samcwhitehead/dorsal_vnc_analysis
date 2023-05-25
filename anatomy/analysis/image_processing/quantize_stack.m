%==========================================================================
% function to use imquantize on a confocal stack and visualize it
%==========================================================================
function segMat = quantize_stack(channelMat, N_levels) 

% perform the multithresh 
levels = multithresh(channelMat, N_levels) ; 
segMat = imquantize(channelMat, levels) ; 

% create a stack of rgb images to display
[nx, ny, nz] = size(channelMat) ; 
rgbMat = uint8(zeros(nx, ny, 3, nz)) ; 

for i = 1:nz
   rgbMat(:,:,:,i) = label2rgb(squeeze(segMat(:,:,i))) ;  
end

implay(rgbMat)

end