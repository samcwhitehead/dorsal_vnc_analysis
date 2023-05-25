%==========================================================================
% copied with minor alterations from Robie, et al 2017 "Mapping the Neural
% Substrates of Behavior." There, they used a lasso regression to fit for
% coefficients that should give thresholds based on a number of image stats
%==========================================================================

function [imnorm,fg_thresh_fit,bg_thresh_fit,imnorm_small,x] = ...
    normalize_channel_flybowl(imgreen,mask, params)

[x,~,imgreen_small] = my_norm_channel_featureVec(imgreen,mask, params);
if isfield(params,'doquad') && params.doquad,
  x = [x;x.^2];
end
fg_thresh_fit = params.fg_coeffs'*x + params.fg_intercept;
bg_thresh_fit = params.bg_coeffs'*x + params.bg_intercept;

fg_thresh_fit = max(params.min_fg_thresh,min(params.max_fg_thresh,fg_thresh_fit));
bg_thresh_fit = max(params.min_bg_thresh,min(params.max_bg_thresh,bg_thresh_fit));
fg_thresh_fit = max(fg_thresh_fit,bg_thresh_fit);

imnorm = min(1,max(0,(imgreen-bg_thresh_fit)/(fg_thresh_fit-bg_thresh_fit)));
imnorm_small = min(1,max(0,(imgreen_small-bg_thresh_fit)/(fg_thresh_fit-bg_thresh_fit)));

end
%--------------------------------------------------------------------------
% my version of the function to get feature vector
function [x,imgreen,imgreen_small] = ...
    my_norm_channel_featureVec(imgreen, mask, params)

% dilate and downsample mask 
mask_small = downsample_stack(mask) ; 
mask_dilated = imdilate(mask_small, params.se_dilate) ; 

% scale the image (box filter + subsample)
%imgreen_smooth = convn(imgreen,params.boxfil,'valid');
imgreen_smooth = imgaussfilt3(imgreen) ; 
imgreen_small = downsample_stack(imgreen_smooth);
[ny,nx,nz] = size(imgreen_small);

% find parts of the image that are padded 0s from registration
is0 = imgreen_small==0;
ccs = bwconncomp(is0);
doignore = false(size(imgreen_small));
for j = 1:numel(ccs.PixelIdxList)
  [iy,ix,iz] = ind2sub([ny,nx,nz],ccs.PixelIdxList{j});
  if any(iy==1) || any(iy==ny) || ...
      any(ix==1) || any(ix==nx) || ...
      any(iz==1) || any(iz==nz)
    doignore(ccs.PixelIdxList{j}) = true;
  end
end
  
mask_bg = 1 - mask_dilated;
mask_bg(doignore) = 0;
mask_fg = double(mask_small);
mask_fg(doignore) = 0;
    
% mean and standard deviation
[mean_bg,var_bg] = weighted_mean_cov(imgreen_small(:),mask_bg(:));
std_bg = sqrt(var_bg);
[mean_fg,var_fg] = weighted_mean_cov(imgreen_small(:),mask_fg(:));
std_fg = sqrt(var_fg);
if params.usegrad 
  % find gradient
  imgreen_mask = imgreen_small;
  imgreen_mask(mask_fg<.5) = 0;
  imgreen_grad = cat(1,zeros(1,nx,nz),diff(imgreen_mask,1,1).^2) + ...
    cat(2,zeros(ny,1,nz),diff(imgreen_mask,1,2).^2);
  isdiff = imgreen_grad > std_bg;
  mean_fg_grad = mean(imgreen_small(isdiff));
  std_fg_grad = std(imgreen_small(isdiff),1);
end
  
% histograms
counts_bg = hist(imgreen_small(:),params.bin_centers_bg,'weights',mask_bg(:));
frac_bg = counts_bg / sum(counts_bg(:));
counts_fg = hist(imgreen_small(:),params.bin_centers_bg,'weights',mask_fg(:));
frac_fg = counts_fg / sum(counts_fg(:));
if params.usegrad,
  counts_fg_grad = hist(imgreen_small(isdiff),params.bin_centers_fg,'weights',mask_fg(:));
  frac_fg_grad = counts_fg_grad / max(1,sum(counts_fg_grad(:)));
end
  
% percentiles
prctiles_bg = weighted_prctile(imgreen_small(:),params.prctiles_compute,mask_bg(:));
prctiles_fg = weighted_prctile(imgreen_small(:),params.prctiles_compute,mask_fg(:));
if params.usegrad
  prctiles_fg_grad = prctile(imgreen_small(isdiff),params.prctiles_compute);
end
  
% construct feature vector
if params.usegrad
  x = [mean_bg;std_bg;mean_fg;std_fg;mean_fg_grad;std_fg_grad;frac_fg(:);frac_bg(:);frac_fg_grad(:);prctiles_fg(:);prctiles_bg(:);prctiles_fg_grad(:)];
else
  x = [mean_bg;std_bg;mean_fg;std_fg;frac_fg(:);frac_bg(:);prctiles_fg(:);prctiles_bg(:)];
end

end

