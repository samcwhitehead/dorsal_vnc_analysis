%==========================================================================
% downsample confocal stack to make plotting/computing less memory
% intensive
%==========================================================================
function imMat_ds = downsample_stack(imMat, dx, dy, dz) 
% make some defaults part of the function for ease of use
if ~exist('dx','var')
    dx = 5 ; 
end
if ~exist('dy','var')
    dy = 5 ; 
end
if ~exist('dz','var')
    dz = 3 ; 
end

% perform downsampling
if size(size(imMat),2) == 4
    imMat_ds = imMat(1:dx:end,1:dy:end,:,1:dz:end) ; 
elseif size(size(imMat),2) == 3
    imMat_ds = imMat(1:dx:end,1:dy:end,1:dz:end) ; 
else
    disp('incorrect image dimensions')
    keyboard
end

end