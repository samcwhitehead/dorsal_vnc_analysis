%==========================================================================
% function to individually normalize channel data from confocal stack
% "lambda_low" and "lambda_high" are the upper and lower thresholds 
%
% *NB: let's assume for now we're working with the downsampled images
%==========================================================================
function channelNorm = normalize_channel(channelMat, mask, debugFlag) 
% try to visualize effects of normalization?
if ~exist('debugFlag','var')
    debugFlag = false ; 
end
% set some parameters
prc_lvl = 99.9 ; % percentile level to use for lower and upper bounds
N_im_lvls = 4 ; % how many intensity levels to decompose image into
gauss_sigma = 1 ; % size for gaussian filter
se = ones(7,7,5) ; % structure element for morphological operations
min_obj_vol = 1000 ; 
defineStacksConstants 

% convert to single array if necessary
if isinteger(channelMat)
    channelMat = single(channelMat).*uint16_to_single ; 
end

% perform top hat transformation to get rid of bg haze
channelMat = imtophat(channelMat, se) ; 

% try to remove some noise
channelMat_smooth = medfilt3(channelMat) ; 
channelMat_smooth = imgaussfilt3(channelMat_smooth,gauss_sigma) ; 

% use downsampled images to make computation easier
channelMat_small = downsample_stack(channelMat_smooth) ;
mask_small = downsample_stack(mask) ; 

% get foreground and background pixels 
foreground_pix = channelMat_small.*mask_small ; 
background_pix = channelMat_small.*~mask_small ; 

max_bg_pix = prctile(background_pix(:), prc_lvl) ; 
%max_fg_pix = prctile(foreground_pix(:), prc_lvl) ; 
max_fg_pix = max(foreground_pix(:)) ; 

%apply normalization
lambda_low = max_bg_pix;
lambda_high = max_fg_pix ; 

%channelNorm = imadjustn(channelMat,[lambda_low, lambda_high],[0,1]) ; 
channelNorm = min(1,max(0,(channelMat-lambda_low)/(lambda_high-lambda_low)));

%channelNorm = imtophat(channelNorm, se) ; 
levels = multithresh(channelNorm,N_im_lvls) ; 
mask_new = (channelNorm > levels(end-1)) ; 
mask_new = bwmorph3(mask_new,'clean') ; 
mask_new = imclose(mask_new,se) ; 
mask_new = imdilate(mask_new, strel('sphere',10)) ; 
% CC = bwconncomp(mask_new,6) ; 
% connVols = cellfun(@(x) length(x), CC.PixelIdxList) ; 

mask_new = bwareaopen(mask_new, min_obj_vol) ; 
%mask_new = imopen(mask_new, se) ; 
channelNormMasked = channelNorm .* mask_new ; 


%test = quantize_stack(channelNorm,2) ; 
% compare mips
if debugFlag
    h_raw = plot_channel_mips(channelMat) ; 
    title('Raw')
    h_norm = plot_channel_mips(channelNorm) ; 
    title('Normalized')
%     h_tophat = plot_channel_mips(channelNorm_tophat) ; 
%     title('Normalized, top hat')
    h_norm_masked = plot_channel_mips(channelNormMasked) ; 
    title('Normalized, Masked')
end
end