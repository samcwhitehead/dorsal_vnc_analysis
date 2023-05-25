% -------------------------------------------------------------------------
% quick script to try to automate the process of registering the segmented
% VNC image from Hiro (label matrix form)
% -------------------------------------------------------------------------
%% define path info
% load seg image
rootPath = 'D:\Fly Imaging\Erica Interneuron Stacks\VNC neuropil\' ;
savePath = fullfile(rootPath, 'pre_process') ;
segPath = fullfile(rootPath, 'segmented_images') ;
confocalFn = 'confocal.tif' ;

overWriteFlag = false ; 

% save name
saveName = 'neuropil_sam_preprocess' ;

% channels per image
N_channels = 1 ;

% target image dimensions
targetDim = [1119, 573, 219] ;

% % load tforms
% tformPathFull = fullfile(rootPath, 'segmented_images', 'old', ...
%     'registered', 'tform_cell.mat') ;
% tform_cell = importdata(tformPathFull) ;

% ---------------------------------------------------------
%% get directory struct for all images
confocalDir = dir(fullfile(rootPath,confocalFn)) ;
segDir = dir(fullfile(segPath, '*.tif')) ;

% combine data directories
dataDir = [confocalDir ; segDir] ;

% ---------------------------------------------------------
%% TIF parameters
tagstruct = struct() ;
tagstruct.ImageLength = targetDim(1) ;
tagstruct.ImageWidth = targetDim(2) ;
tagstruct.Photometric = 1 ; % 1 = "min is black"
tagstruct.BitsPerSample = 8; % saving as 8bit image
tagstruct.SamplesPerPixel = 1 ; % aka number of channels
tagstruct.Compression = 1 ; %1 = "none"
tagstruct.PlanarConfiguration = 1 ; % "Chunky" configuration. who fucking knows?
tagstruct.Software = 'MATLAB';

% ---------------------------------------------------------
%%  constants for transformation
% rotation
M_rot = [0.4002, -0.9164, 0, 0 ; 0.9164, 0.4002, 0, 0 ; 0, 0, 1.0,  0 ; ...
    0, 0, 0, 1.0] ;
tform_rot = affine3d(M_rot) ;

% removing empty slices
good_z_ind = 62:182 ;  % (these are slices with VNC in them)

% -------------------------
% crop info
xScale = 1.0/1.5 ; % NB: scaling determined from point cloud registration (but tweaked a lil)
yScale = 1.0/1.5 ;

cropHeight = yScale*targetDim(1) ;
cropWidth = xScale*targetDim(2) ;
cropDepth = 120 ; % NB: this is the image size after cropping z stacks
% center of mass in rotated image
%{
segBW = imbinarize(segIm) ;
rprops = regionprops3(segBW, 'Centroid', 'Volume') ;
[~, max_vol_ind] = max(rprops.Volume) ;
cm_seg = rprops.Centroid(max_vol_ind,:) ;
%}
cm_seg = [815.5243, 836.3878 + 20.0 , 65.9078] ;

% use cm and cropWidth/Height to define box for cropping
xmin = round(cm_seg(1) - cropWidth/2) ;
ymin = round(cm_seg(2) - cropHeight/2) ;
zmin = 1 ;

cuboid = [xmin, ymin, zmin, cropWidth, cropHeight, cropDepth] ;

% ---------------------------------------------------
%% loop over images
N_im = length(dataDir) ;

for ind = 1:N_im
    % check if processed already exists
    savePathFull = fullfile(savePath, ...
        [saveName '_' num2str(ind,'%02d') '.tif']) ;
    
    if exist(savePathFull,'file') && ~overWriteFlag
       fprintf('Already processed %s -- skipping \n', savePathFull)
       continue
    end
    % load data
    dataPathFull = fullfile(dataDir(ind).folder, dataDir(ind).name) ;
    segIm = myReadTiff(dataPathFull, N_channels) ;
    
    % -------------------------------------------
    % first transform: apply rotation about z axis
    segIm = imwarp(segIm, tform_rot) ;
    
    % -------------------------------------------
    % remove slices with nothing in them
    segIm = segIm(:,:,good_z_ind) ;
    
    % --------------------------------------------
    % crop area around image in xy
    segIm = imcrop3(segIm, cuboid) ;
    
    % ------------------------------------------------
    % resize image to the proper dimensions
    % NB: extrapolation method depends on whether we're dealing with true
    % grayscale image (original confocal data) or basically binary/label
    % matrix (individual neuropil images)
    if ind == 1
        segIm = imresize3(segIm, targetDim) ;
    else
        segIm = imresize3(segIm, targetDim, 'nearest') ;
    end
    % -------------------------------------------------
    % convert image to 8 bit
    segIm = im2uint8(segIm) ;
    
    % --------------------------------------------------------------
    % adjust image pixel intensity (only on neuropil images, which are
    % essentially binary anyway)
    if ind ~= 1
        segIm = imadjustn(segIm) ;
    end
    
    % -------------------------------------------------
    %% save image
    
    % initialize Tiff object
    t = Tiff(savePathFull, 'w') ;
    
    % loop over images and write each
    imDepth = size(segIm,3) ;
    for k = 1:imDepth
        % need to set tags for each directory?
        setTag(t, tagstruct) ;
        
        % write image
        write(t, squeeze(segIm(:,:,k))) ;
        
        % advance to next slice (unless we're at end)
        if k ~= imDepth
            writeDirectory(t) ;
        end
    end
    
    % close tif
    close(t)
    
    fprintf('Completed pre-processing %d/%d images \n', ind, N_im)
end