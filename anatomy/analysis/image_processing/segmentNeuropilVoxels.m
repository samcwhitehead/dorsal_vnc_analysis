% -------------------------------------------------------------------------
% script to take tif image with segmented neuropil regions in confocal
% stack (from Hiro) and save separate images with each region
% -------------------------------------------------------------------------
%% path info
dataPath = 'D:\Fly Imaging\Erica Interneuron Stacks\VNC neuropil\' ;

% flags and misc. params
saveFlag = true ;        % save results?
overWriteFlag = false ;  % over write previous files
% dilateFlag = true ;     % make masks a little larger to fill in some of the gaps?
imFileExt = '.tif' ;     % what file types are images?
debugBWFlag = false ;

% where to look for data
dataPath = fullfile(dataPath, 'segmented_images', 'registered') ;

% save info
savePath = dataPath ;
if ~exist(savePath,'dir')
    mkdir(savePath)
end

% -----------------------
% params
imDim = [1119, 573, 333] ; 
shrinkFactor = 1 ; % amount to shrink volume by to speed up computation
%voxSize = [0.61, 0.61, 0.731] ; %note: this is just an estimate--don't have true dimensions
cleanFlag = true ; % clean BW image?
debugBoundaryFlag = false ; % check output?
N_channels = 1 ; % number of channels in image
smoothFactor = 0.1 ;
boundaryConv = 0.9 ; % makes tight boundary
minObjSize = 5e3 ;

% structural element
se = strel('sphere',5) ;

% corresponding neuropil labels for each pixel value
neuropilLabels = {'LTct',... % 1
    'WTct', ... %2
    'AMNp', ... %3
    'mVAC', ... %4
    'VAC1', ... %5
    'NTct', ... %6
    'T1LegNp', ... %7
    'T2LegNp', ... %8
    'T3LegNp', ... %9
    'ANm', ... %10
    'HTct', ... %11
    'VAC2', ... %12
    'VAC3', ... %13
    'IntTct',... %14
    'IntTct1',...
    'IntTct2',...
    'IntTct3',...
    'LTct1',...
    'LTct2'} ; 

% ---------------------------------------------------
%% get directory for neuropil images
neuropilDir = dir(fullfile(dataPath, ['*' imFileExt])) ;
N_np = length(neuropilDir) ;

% ---------------------------------------------------
%% loop through neuropil images
for i = 1:N_np
    % ----------------------------------------------------------------
    % get path to current neuropil image
    fn = neuropilDir(i).name ;
    path = neuropilDir(i).folder ;
    [~, regionLabel, ~] = fileparts(fn) ; % get name of neuropil from filename
    
    % get save path and check if we've already done this
    savePathFull = fullfile(savePath, regionLabel) ;
    if exist([savePathFull '.mat'],'file') && ~overWriteFlag
        fprintf('Already processed %s -- continuing \n', regionLabel)
        continue
    end
    
    % --------------------------------------------------------------
    %% load image of current neuropil region
    fprintf('Loading image for %s ... \n', regionLabel)
    imMat = myReadTiff(fullfile(path, fn), N_channels) ;
    
    % check to see if image is "VNC" -- this is actually the original
    % confocal stack, so require different processing
    if strcmpi(regionLabel, 'VNC')
        filt_sigma = 1.5 ;
    else
        filt_sigma = 1.0 ;
    end
    
    % ------------------------------------------------------
    %% process image
    fprintf('Processing image for %s ... \n', regionLabel)
    
    % resize image so that voxels have equal dimensions (i.e. are cubes)
    imMat = imresize3(imMat, imDim, 'nearest') ; 
    
    % filter image to account for some of the weirdness arising from reg.
    imMatFilt = imgaussfilt3(imMat, filt_sigma) ;  
    
    % get multiple threshold values -- use lower one
    thresh_levels = multithresh(imMatFilt,2) ; 
    
    % binarize image using lower thresh
    regionBW = (imMatFilt >= thresh_levels(1)) ; % binarize
    
    % remove noisy junk
    if strcmpi(regionLabel, 'VNC')
        % VNC should only have one object
        regionBW = imopen(regionBW, se) ; 
        CC = bwconncomp(regionBW) ;
        cc_sizes = cellfun(@(y) length(y), CC.PixelIdxList) ; 
        max_cc_size = max(cc_sizes) ; 
        
        regionBW = bwareaopen(regionBW, max_cc_size - 1) ; 
        
    else
        % otherwise just remove small globs
        regionBW = bwareaopen(regionBW, minObjSize) ;
    end

    % morphological closing
    regionBW = imclose(regionBW,se) ;
    
    % extra VNC processing
    if strcmpi(regionLabel, 'VNC')
        % dilate a little so the VNC nicely encloses neuropil regions
        regionBW = imdilate(regionBW,strel('sphere',1)) ;
    end
    
    % if we're dealing with the abdominal neuromere, need to get convex
    % hull (the registration seems to chunk up the image in a weird way)
    if strcmpi(regionLabel, 'ANm')
        % get convex hull voxel image and bounding box (the convex image
        % output is cropped to bounding box region)
       rprops = regionprops3(regionBW,'BoundingBox', 'ConvexImage') ; 
       bbox = rprops.BoundingBox ; 
       ch = rprops.ConvexImage{1} ; 
       
       % insert convex hull into image
       indx = round(bbox(1)):(round(bbox(1) + bbox(4)) - 1) ; 
       indy = round(bbox(2)):(round(bbox(2) + bbox(5)) - 1) ; 
       indz = round(bbox(3)):(round(bbox(3) + bbox(6)) - 1) ; 
       [xx, yy, zz] = meshgrid(indx, indy, indz) ; 
       
       ch_idx = sub2ind(size(regionBW), yy(:), xx(:), zz(:)) ; 
       regionBW(ch_idx) = ch(:) ; 
       
       % keyboard
    end
      
    % check that binarization went okay?
    if debugBWFlag
        figure ;
        imshow3D(regionBW)
        
        figure
        imshowpair(max(imMat,[],3), max(regionBW,[],3))
        keyboard
        
        close all
    end
    % -----------------------------------------------------
    %% find connected components in image
    fprintf('Converting image of %s to matlab data... \n', regionLabel)
    CC = bwconncomp(regionBW) ;
    numObjects = CC.NumObjects ;
    

    % initialize output storage
    vox_coord_cell = cell(numObjects,1) ;  % boundary coordinates
    vox_boundary_cell = cell(numObjects,1) ; % boundary triangulation
    bw_cell = cell(numObjects,1) ; % binary image (full)
    coord_ind_cell = cell(numObjects,1) ; % coordinates for full BW image, in index form
    
    % loop through objects and get coordinates for each
    for j = 1:numObjects
        % get image for CC
        regionBW_curr = false(size(regionBW)) ;
        regionBW_curr(CC.PixelIdxList{j}) = true ;
        
        % get boundary points and mesh
        [vox_coords, vox_boundary] = getVoxelBoundaryPoints(regionBW_curr, ...
            smoothFactor, boundaryConv, cleanFlag, debugBoundaryFlag) ;
        vox_coord_cell{j} = vox_coords ;
        vox_boundary_cell{j} = vox_boundary ;
        bw_cell{j} = regionBW_curr ;
        coord_ind_cell{j} = find(regionBW_curr(:)) ;
        
    end
    
    % ---------------------------------------------------
    %% save?
    if saveFlag
        save(savePathFull,'vox_coord_cell','vox_boundary_cell','bw_cell')
        save([savePathFull '_coords.mat'], 'coord_ind_cell')
    end
    
    fprintf('Completed %s (%d/%d) \n',regionLabel, i, N_np)
end

% ----------------------------------------------------------------
%% final thing: combine VAC1-3 into one VAC file
% (if it doesn;t already exist)
vacSavePath = fullfile(savePath, 'VAC.mat') ;
vacCoordsSavePath = fullfile(savePath, 'VAC_coords.mat') ;

% if VAC file already exists, do nothing. otherwise combine files
if ~exist(vacSavePath,'file') || ~exist(vacCoordsSavePath,'file')
    fprintf('Combining VAC files ... \n')
    % get directory for VAC1-3 data
    vacDir = dir(fullfile(dataPath, 'VAC*.mat')) ; 
    regionLabel = 'VAC' ; 
    savePathFull = fullfile(savePath, regionLabel) ;
    
    % deal with coords and image data separately
    coords_idx = arrayfun(@(x) contains(x.name, '_coords'), vacDir) ; 
    
    % --------------------------------------------------------
    % COORDS
    vacCoordsDir = vacDir(coords_idx) ; 
    for q = 1:length(vacCoordsDir)
        load_data = importdata(fullfile(vacCoordsDir(q).folder, ...
            vacCoordsDir(q).name)) ;
        if q == 1
            coord_ind_cell = load_data ;
        else
            coord_ind_cell = [coord_ind_cell ; load_data];
        end
    end
    
    % save coords?
    if saveFlag
        save([savePathFull '_coords.mat'], 'coord_ind_cell')
    end
    
    % ------------------------------------------------------
    % IMAGE AND BOUNDARY
    vacImDir = vacDir(~coords_idx) ; 
    for q = 1:length(vacImDir)
        load_data = importdata(fullfile(vacImDir(q).folder, ...
            vacImDir(q).name)) ;
        if q == 1
            vox_coord_cell = load_data.vox_coord_cell ; 
            vox_boundary_cell = load_data.vox_boundary_cell ; 
            bw_cell = load_data.bw_cell ; 
        else
            vox_coord_cell = [vox_coord_cell ; load_data.vox_coord_cell] ; 
            vox_boundary_cell = [vox_boundary_cell ; ...
                load_data.vox_boundary_cell] ; 
            bw_cell = [bw_cell; load_data.bw_cell] ; 
        end
    end
    
    % save image and boundary
    if saveFlag
        save(savePathFull,'vox_coord_cell','vox_boundary_cell','bw_cell')
    end
end


