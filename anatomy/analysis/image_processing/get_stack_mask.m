%==========================================================================
% create binary mask for a tif stack that gives the full volume of the CNS
% region. NB: this assumes that the brightfield info is stored in channel 1
%==========================================================================
function mask = get_stack_mask(imMat, debugFlag) 

% debug results?
if ~exist('debugFlag','var')
    debugFlag = false ; 
end
% structuring element for morpho operations
se_close = ones(7,7,5) ; 
se_open = ones(5,5,3) ; 
% define some params
defineStacksConstants

% read out channel 1 (assumed to be brightfield) and convert to normalized
% double
bfield_stack = squeeze(imMat(:,:,BFIELD,:)) ; 
bfield_stack = double(bfield_stack)./double(max(bfield_stack(:))) ; 

% apply gaussian filter to smooth out some of the noise
bfield_stack_filt = imgaussfilt(bfield_stack) ; 

% binarize the image using otsu's method
level = graythresh(bfield_stack_filt) ;
bfield_stack_bw = imbinarize(bfield_stack_filt, level) ; 

% now perform morphological operations to clean the image up
mask = bwmorph3(bfield_stack_bw, 'clean') ; 
% close image to fill gaps
mask = imclose(mask,se_close) ;
% fill interior voxels
mask = bwmorph3(mask, 'fill') ; 
% image opening?
mask = imopen(mask,se_open) ; 

if debugFlag
    implay(mask)
end


end