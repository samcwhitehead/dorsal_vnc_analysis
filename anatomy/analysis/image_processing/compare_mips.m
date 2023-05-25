function [score, image_comp] = compare_mips(im_target, im_search, cmap)

plotFlag = true ;
%filtFlag = true ; 
KernelSigma = 0.01 ; % in units of double precision image 

% convert rgb to indices (which represent height
im_target_ind = rgb2ind(im_target,cmap) ;
im_search_ind = rgb2ind(im_search,cmap) ;

%im_target_ind = imcomplement(im_target_ind) ;
%im_search_ind = imcomplement(im_search_ind) ;

% covnert to double precision
im_target_double = im2double(im_target_ind) ;
im_search_double = im2double(im_search_ind) ;

% black is not in the mip colormap, so needs to be removed
bg_pix = mode(im_target_double(:)) ;
im_target_double(im_target_double == bg_pix) = nan ;
im_search_double(im_search_double == bg_pix) = nan ;


if (0)
    % just to check color -> z conversion
    figure ;
    
    subplot(1,2,1)
    imshow(im_target_double)
    colormap(gca,hot)
    
    subplot(1,2,2)
    imshow(im_search_double)
    colormap(gca,hot)
end

% calculate gaussian kernel distance 
im_dist = double(im_target_double - im_search_double) ;
im_sq_dist = im_dist.^2 ;
im_gauss_dist = exp(-1*(1/2*KernelSigma^2).*im_sq_dist) ;

if plotFlag
    figure ;
    subplot(1,2,1)
    imshow(im_gauss_dist,[])
    colormap(gca,hot)
    
    subplot(1,2,2)
    imshowpair(im_target_double,im_search_double)
end

score = nan ; 
image_comp = im_gauss_dist ;

end