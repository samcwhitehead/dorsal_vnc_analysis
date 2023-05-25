% -------------------------------------------------------------------------
% attempt at processing image stacks using canned MATLAB functions
% -------------------------------------------------------------------------
function BW = binarize_confocal_stack(I, se_rad, ...
    removeASFlag, plotFlag)
% --------------------
%% params and inputs
if ~exist('se_rad','var') || isempty(se_rad)
    se_rad = 3 ; % radius for spherical structuring element used to process
end
if ~exist('removeASFlag','var') || isempty(removeASFlag)
    removeASFlag = false ; % remove abdominal segment or no? (bottom of posterior-most part of VNC)
end
if ~exist('plotFlag','var') || isempty(plotFlag)
    plotFlag = false ; % show results?
end

N_levels = 3 ; % how many levels to quantize image into
cc_z_thresh = 1 ; % threshold for connected component object size (in z score)
se = strel('sphere',se_rad) ;
high_prctile = 98 ;

% ------------------------------------------------------------
%% get MIP of image and use to determine values for imadjust
% quantize image into N levels
MIP = max(I, [], 3) ;
MIP = imgaussfilt(MIP) ; 
[levels, thresh_metric] = multithresh(MIP, N_levels) ;

% get low_in and high_in for imadjust. for low, use multithresh levels. for
% high use percentile of MIP image
if isa(MIP,'integer')
    low_in = double(levels(1))/double(intmax(class(MIP))) ;
    high_in = double(prctile(MIP(:), high_prctile))/double(intmax(class(MIP))) ;
else
    low_in = double(levels(1)) ;
    high_in = double(prctile(MIP(:), high_prctile)) ;
end

% make sure that the limits are in increasing order
while (high_in < low_in)
    high_prctile = high_prctile + 0.75  ;
    if isa(MIP,'integer')
        high_in = double(prctile(MIP(:), high_prctile))/double(intmax(class(MIP))) ;
    else
        high_in = double(prctile(MIP(:), high_prctile)) ;
    end
end
%  adjust stack image contrast
J = imadjustn(I, [low_in, high_in], [0, 1]) ;
J = imgaussfilt3(J) ; 
% ------------------------------------------------------------
%% binarize and clean up image stack
%BW0 = imbinarize(J) ;
%BW0 = (J > levels(end)) ;
[levels_3d, thresh_metric_3d] = multithresh(J, N_levels) ; 
if (thresh_metric_3d == 0)
    BW0 = imbinarize(J) ; 
else
    BW0 = (J > levels_3d(end)) ; 
end

% close image to fill in gaps
BW1 = imclose(BW0, se) ;

% get distribution of connected components to remove small globs
CC = bwconncomp(BW1) ;
objSizes = cellfun(@(y) length(y), CC.PixelIdxList) ;

% remove objects whose z scored size is smaller thanthreshold
% (set as 'cc_z_thresh' above)
minObjSize = round(cc_z_thresh*std(objSizes)) ;
BW2 = bwareaopen(BW1, minObjSize) ;

% fill any remaining wholes in image
BW = bwmorph3(BW2, 'fill') ;

% ----------------------------------------------------------------
%% remove abdominal segment?
if removeASFlag
    imHeight = size(BW,1) ;
    ind1 = round(0.8*imHeight) ;
    ind2 = imHeight ;
    BW(ind1:ind2,:,:) = false ;
end

%% ----------------------------------------------------------
%% show results?
if plotFlag
    % original
    volumeViewer(I) ;
    % output
    volumeViewer(BW) ;
end

end