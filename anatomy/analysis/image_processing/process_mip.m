%--------------------------------------------------------------------------
% take a mip as input and crop/thresh/mask it for comparison
%--------------------------------------------------------------------------
function [mip_im_masked, mip_mask, mip_im] = ...
    process_mip(driver,rootPath,plotFlag1,plotFlag2)

if nargin < 2
    rootPath = 'C:\Users\samcw\Dropbox\mips\Erica\MIP\' ; % where to find mips
    plotFlag1 = false ; % show montage of all the processing stuff
    plotFlag2 = false ; % show result of masking
elseif nargin < 3
    plotFlag1 = false ; % show montage of all the processing stuff
    plotFlag2 = false ; % show result of masking
elseif nargin < 4
    plotFlag2 = false ; % show result of masking
end

% general mip info
%driver = 'SS37252' ;
fileInfoStr = 'F_02_warp' ;
headerLength = 90 ;

% structural element for morphological operations
morphSize = 10 ;
se = strel('disk',morphSize) ;
conn = 4 ; % required connectivity
CC_thresh = 200 ; % minimum size of connected components
sigma = 1 ; 

multiThreshFlag = false ; 
N_thresh = 3 ; 

%figPosition = 1.0e+03*[0.0817    0.2357    1.1740    0.3020] ;
im2plot_cell = cell(1) ;
title_cell = cell(1) ;
cc = 1 ;
%--------------------------------------------------------------------------
%% load image
dataDir = dir([rootPath '*' fileInfoStr '*.tif']) ;
driver_ind = arrayfun(@(x) contains(x.name,driver(2:end)), dataDir) ;
if sum(driver_ind) ~= 1
    disp('failed to find single image')
    keyboard
end
dataPath = fullfile(dataDir(driver_ind).folder, dataDir(driver_ind).name) ;
mip_im = importdata(dataPath) ;

%--------------------------------------------------------------------------
%% crop image to remove header
mip_im = mip_im((headerLength+1):end, :, :) ;

im2plot_cell{cc} = mip_im ;
title_cell{cc} = 'cropped im' ;
cc = cc + 1 ;
%--------------------------------------------------------------------------
%% convert to hsv to threshold without affecting depth information
mip_im_hsv = rgb2hsv(mip_im) ;

im2plot_cell{cc} = mip_im_hsv ;
title_cell{cc} = 'hsv conversion' ;
cc = cc + 1 ;


% just get value--shows neural expression well
mip_im_value = squeeze(mip_im_hsv(:,:,3)) ;

im2plot_cell{cc} = mip_im_value ;
title_cell{cc} = 'value of im' ;
cc = cc + 1 ;

%--------------------------------------------------------------------------
%%  gaussian filter of value image
mip_im_val_filt = imgaussfilt(mip_im_value,sigma) ; 
im2plot_cell{cc} = mip_im_val_filt ;
title_cell{cc} = 'gaussian filter' ;
cc = cc + 1 ;

%--------------------------------------------------------------------------
%%  binarize image using otsu method
if multiThreshFlag
    thresh = multithresh(mip_im_val_filt,N_thresh);
    mip_bw = imbinarize(mip_im_val_filt, thresh(2)) ;
%     seg_I = imquantize(mip_im_val_filt,thresh);
%     RGB = label2rgb(seg_I);
%     figure;
%     imshowpair(mip_im_value,RGB,'montage')
%     axis off
%     title('RGB Segmented Image')
else
    level = graythresh(mip_im_val_filt) ;
    mip_bw = imbinarize(mip_im_value, level) ;
    mip_filt_bw = imbinarize(mip_im_val_filt, level) ;
    %figure ; imshowpair(mip_bw, mip_filt_bw,'montage')
end


im2plot_cell{cc} = mip_bw ;
title_cell{cc} = 'binarized w otsu thresh' ;
cc = cc + 1 ;

%--------------------------------------------------------------------------
%% morphological operations
%mip_bw1 = medfilt2(mip_bw) ;
mip_bw2 = imclose(mip_bw,se) ;

im2plot_cell{cc} = mip_bw2 ;
title_cell{cc} = 'morph. bw' ;
cc = cc + 1 ;

%--------------------------------------------------------------------------
%% connected components
CC = bwconncomp(mip_bw2,conn) ;
CC_sizes = cellfun(@(x) length(x), CC.PixelIdxList) ;
good_CC_idx = find(CC_sizes >= CC_thresh) ;

mip_mask = false(size(mip_bw2)) ;
for k = 1:length(good_CC_idx)
    mip_mask(CC.PixelIdxList{good_CC_idx(k)}) = true ;
end
%mip_bw3 = imfill(mip_bw3,'holes') ;

im2plot_cell{cc} = mip_mask ;
title_cell{cc} = 'conn comp' ;
cc = cc + 1 ;

%% plot processing results?
if plotFlag1
    h_montage = figure('PaperPositionMode','auto') ; %,'Position',figPosition);
    for i = 1:(cc-1)
        subplot(1,(cc-1),i)
        imshow(im2plot_cell{i})
        title(title_cell{i})
    end
end

%% use connected component as mask
mip_im_masked = mip_im ;
for j = 1:3
    mip_im_channel = mip_im(:,:,j) ;
    mip_im_channel = mip_im_channel.*uint8(mip_mask) ;
    mip_im_masked(:,:,j) = mip_im_channel ;
end

if plotFlag2
    figure ;
    subplot(1,3,1)
    imshow(mip_im)
    title('Original Im')
    
    subplot(1,3,2)
    imshow(mip_im_masked)
    title('Masked Im')
    
    subplot(1,3,3)
    imshowpair(mip_im, mip_im_masked)
    title('Overlay')
end
