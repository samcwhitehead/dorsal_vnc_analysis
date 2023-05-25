%--------------------------------------------------------------------------
% convert rgb values of mip to linear z coordinates
%--------------------------------------------------------------------------

% which image to take (shouldn't matter)
driver = 'SS37252' ;
fileInfoStr = 'F_02_warp' ;
headerLength = 90 ;

autoSelectFlag = true ;
saveFlag = true ; 
plotFlag = false ; 

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

%% select area of image with color scale
mip_header = mip_im(1:headerLength,:,:) ;

if autoSelectFlag
    % find areas where the image header is either black or white
    white_pix_val = max(mip_header(:)) ;
    black_ind = find(mip_header(:,:,1) == 0 & mip_header(:,:,2) == 0 ...
        & mip_header(:,:,3) == 0) ;
    white_ind = find(mip_header(:,:,1) == white_pix_val & ...
        mip_header(:,:,2) == white_pix_val ...
        & mip_header(:,:,3) == white_pix_val) ;
    
    % make mask excluding these pixels. This should mostly leave the
    % colorbar. Use connected comps to clean up
    cb_mask = true(size(mip_header,1), size(mip_header,2)) ; 
    cb_mask(black_ind) = false ; 
    cb_mask(white_ind) = false ; 
    
    CC = bwconncomp(cb_mask) ; 
    CC_sizes = cellfun(@(x) length(x), CC.PixelIdxList) ;
    [~, CC_max_ind] = max(CC_sizes) ; 
    
    cb_mask2 = false(size(cb_mask)) ; 
    cb_mask2(CC.PixelIdxList{CC_max_ind}) = true ; 
  
    mask_corners = corner(cb_mask2);
    rect2 = [min(mask_corners(:,1)), min(mask_corners(:,2)), ...
        max(mask_corners(:,1)) - min(mask_corners(:,1)) - 1, ...
        max(mask_corners(:,2)) - min(mask_corners(:,2)) - 1 ] ; 
    if plotFlag
        figure ;
        %subplot(2,1,1)
        imshow(cb_mask2)
        hold on
        plot(mask_corners(:,1),mask_corners(:,2),'r*');
       
    end
else
    rect2 = [256, 2, 255, 31] ;
end

% check that the x range matches an integer size! (the bar is oriented to
% vary along the x axis, and should have 256 values for an 8-bit image,
% etc.
if round(log2(rect2(3)+1))/log2(rect2(3)+1) ~= 1
    disp('under construction, but color range is off!')
    keyboard ; 
end

% crop image
mip_cb = imcrop(mip_header,rect2) ;
if plotFlag
    figure ;
    imshow(mip_cb)  
end

%% define the colormap and save
% take midline to be safe
mid_y = round(rect2(4)) ; 
cmap = squeeze(mip_cb(mid_y, :, :)) ;
cmap = double(cmap)./255 ; 

if saveFlag
    mFileName = mfilename('fullpath') ;
    [savePath,~,~] = fileparts(mFileName) ;
    save(fullfile(savePath, 'MIP_colormap.mat'),'cmap') ;
end


