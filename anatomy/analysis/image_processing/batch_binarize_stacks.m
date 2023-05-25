% -------------------------------------------------------------------------
% script to loop through segmented interneuron stacks from Erica, binarize
% them, resize them, and then save the logical arrays. the hope is that
% this will speed up overlap computations, since much of the overhead seems
% to be i/o
% -------------------------------------------------------------------------
%% data path info
dataRoot = 'D:\Fly Imaging\Erica Interneuron Stacks\VUM\' ;
%dataRoot = 'D:\Fly Imaging\Erica Interneuron Stacks\DN\' ;
tiffDir = dir([dataRoot '*.tif']) ;
nrrdDir = dir([dataRoot '*.nrrd']) ;
imDir = vertcat(nrrdDir, tiffDir) ;
N_im = length(imDir) ;

% save over? plot?
overWriteFlag = false ;
debugFlag = true ;
extraProcessFlag = false ; 

dataType = '' ; % 'DN' | 'DN_mask' | '' 
% should we try to symmetrize?
symFlag = true ;   % will auto-set to false for VUMs


% make sure dataType is set
if isempty(dataType)
   dataRootSplit = strsplit(dataRoot,'\') ; 
   temp_empty_idx = cellfun(@(y) isempty(y), dataRootSplit) ; 
   dataRootSplit = dataRootSplit(~temp_empty_idx) ; 
   dataType = dataRootSplit{end} ; 
end
% ---------------------------------------------------------------
% figure out what kind of data we're dealing with
% dataRootSplit = strsplit(dataRoot,'\') ; 
% empty_cell_idx = cellfun(@(y) isempty(y), dataRootSplit) ; 
% dataRootSplit = dataRootSplit(~empty_cell_idx) ; 

switch dataType
    case 'DN'
        symThresh = 0.05 ; % 0.3 for interneurons, 0.05 for DNs
        
        % try to remove abdominal segment cells? they show up a lot...
        removeASFlag = true ;
        AS_row_start = 905 ;
        
        % image info
        N_channels = 3 ; % the segmented images are just one channel (green)
        green_channel_ind = 2 ;
        DNFlag = true ;  % are we looking at DNs? if so, remove CCs that don't look like they extend to the neck
        DN_top_dist_thresh = 50 ; % lowest we'll allow the top of a CC to be considered a DN
        %NB : the above heuristic could give trouble if a single neuron gets
        %broken up in the processing!
        pixScaleFactor = 1 ; 
        
    case 'DN_mask' 
        % in this case, we're dealing with segmented images
        symThresh = 0.05 ; % 0.3 for interneurons, 0.05 for DNs
        
        % try to remove abdominal segment cells? they show up a lot...
        removeASFlag = false ;
        AS_row_start = [] ;
        
        % image info
        N_channels = 1 ; % the segmented images are just one channel (green)
        green_channel_ind = [] ;
        DNFlag = false ;  % are we looking at DNs? if so, remove CCs that don't look like they extend to the neck
        DN_top_dist_thresh = [] ; 
        
        pixScaleFactor = 255 ; 
    otherwise
        % in this case, we're dealing with segmented images
        symThresh = 0.3 ; % 0.3 for interneurons, 0.05 for DNs
        
        % try to remove abdominal segment cells? they show up a lot...
        removeASFlag = false ;
        AS_row_start = [] ;
        
        % image info
        N_channels = 1 ; % the segmented images are just one channel (green)
        green_channel_ind = [] ;
        DNFlag = false ;  % are we looking at DNs? if so, remove CCs that don't look like they extend to the neck
        DN_top_dist_thresh = [] ; 
        
        pixScaleFactor = 1 ; % 255 
        
        % if we're looking at VUMs, no need to symmetrize
        if strcmpi(dataType,'VUM')
           symFlag = false ;  
        end
end
    
% processing params
maxLevel = 0.075 ;
minObjSize = 1e4 ;
voxDimDefault = [0.461, 0.461, 0.7] ;

% make directory for storing binarize/resized images
bwNonSymPath = fullfile(dataRoot, 'binarized_non_sym') ; 
bwSymPath = fullfile(dataRoot, 'binarized_new') ; 
if symFlag || strcmpi(dataType,'VUM')
    bwPath = bwSymPath ;
else
    bwPath = bwNonSymPath ; 
end
if ~exist(bwPath, 'dir')
    mkdir(bwPath)
end


% indices of stacks where binarization fails
fail_idx = false(N_im, 1) ;
% -------------------------------------------------------------------------
%% get image scale
% since the images are sometimes different sizes, find the smallest
% dimensions. we'll then resize all images to that scale
if exist(fullfile(dataRoot,'imDimensions.mat'),'file') ~= 2
    heightList = zeros(N_im,1) ;
    widthList = zeros(N_im,1) ;
    depthList = zeros(N_im,1) ;
    
    for i = 1:N_im
        % first get size in voxels
        fn_curr = fullfile(imDir(i).folder,imDir(i).name) ;
        imgInfo = imfinfo(fn_curr) ;
        
        heightCurr = imgInfo(1).Height ;
        widthCurr = imgInfo(1).Width ;
        N_images = length(imgInfo) ;
        depthCurr = round(N_images/N_channels) ;
        
        % then use the real-space voxel sizes to get dimensions
        % corresponding to equal space
        voxDimVec = getVoxelDimensions(fn_curr) ;
        if any(isnan(voxDimVec))
            voxDimVec = voxDimDefault ;
        end
        
        minVoxSize = min(voxDimVec) ;
        scaleVec = voxDimVec./minVoxSize ;
        
        heightList(i) = scaleVec(1)*heightCurr ;
        widthList(i) = scaleVec(2)*widthCurr ;
        depthList(i) = scaleVec(3)*depthCurr ;
        
    end
    
    imHeight = min(heightList) ;
    imWidth = min(widthList) ;
    imDepth = min(depthList) ;
else
    load(fullfile(dataRoot,'imDimensions.mat'))
end

% -----------------------------------------------------------
%% loop over images and process them
for ind = 1:N_im
    tic
    
    fprintf('Processing %d / %d stacks \n', ind, N_im)
    % get filename for image
    dataFilename = fullfile(imDir(ind).folder, imDir(ind).name) ;
    [~, fn, ext] = fileparts(dataFilename) ;
    % set save path
    bwSavePath = fullfile(bwPath, [fn '_bw.mat']) ;
    
    % check if we've already converted this one
    if (~overWriteFlag) && exist(bwSavePath,'file')
        bwMat = importdata(bwSavePath) ; 
        if sum(bwMat(:)) > 1
            fprintf('Already analyzed file: %s \n', dataFilename)
            continue
        else
            fprintf('Already analyzed file: %s , but need to redo... \n',...
                dataFilename)
        end
    end
    
    % if we're looking for a symmetrized image, and already have the non
    % sym version, look for that too
    bwNonSymSavePath = fullfile(bwNonSymPath, [fn '_bw.mat']) ;
    if symFlag && exist(bwNonSymSavePath,'file') 
        % load (potentially) unilateral image
       bwMat = importdata( bwNonSymSavePath) ;  
       
       % check if it's actually unilateral
       symCheck = sum(bwMat & flip(bwMat,2),'all')/sum(bwMat(:)) ;
        %disp(symCheck)
        
        % if it seems unilateral, symmetrize
        if symCheck < symThresh
            bwMat = bwMat | flip(bwMat,2) ;
        end
        
        % save symmetrized version (and plot?)
        if debugFlag
            % initialize figure
            figUnits = 'inches' ;
            figPosition = [5.4583, 1.8333, 6.2917, 7.8125] ;
            h_IN = figure('PaperPositionMode','auto','MenuBar','none',...
                'ToolBar','none','DockControls','off','Units',figUnits,...
                'OuterPosition',figPosition) ;
            
            % initialize axis
            ax = gca ;
            hold(ax,'on')
            
            % create a parent hgtform
            parent = hgtransform('Parent',ax);
            
            % plot params
            vncAlpha = 0.05 ;
            neuronAlpha = 1.0 ;
            inputColor = lines(1) ;
            imSize = size(bwMat) ;
            
            % draw VNC background and neuron
            drawVNCBackgroundToAxis(ax, h_IN, [], [], vncAlpha, parent) ;
            ax = drawNeuronToAxis(ax, bwMat, inputColor, neuronAlpha, ...
                imSize, parent) ;
            
            exportgraphics(h_IN, fullfile(bwPath, [fn '_bw_plot.png']),...
                'Resolution', 500)
            close all
        end
        save(bwSavePath,'bwMat')
        
        fprintf('Using unilateral image to make symmetric image ... \n')
        toc
        continue
    end
    % -----------------------------------------------------------
    %% load image
    if strcmp(ext,'.tif')
        imMat  = myReadTiff(dataFilename, N_channels) ;
    elseif strcmp(ext,'.nrrd')
        [imMat, ~] = nrrdread(dataFilename) ;
    else
        fprintf('Invalid file extension: %s \n', ext)
        fail_idx(ind) = true ;
        continue
    end
    
    % ---------------------------------------------------------------
    %% if it's a multi-channel image, grab only channel we care about
    % (usually green)
    if (N_channels > 1)
        % if it's multichannel, then probably not segmented, and we need to
        % work to get binary image out
        imMat = squeeze(imMat(:,:,green_channel_ind,:)) ;
        
        % still need to resize image 
        imMat = imresize3(imMat,[imHeight, imWidth, imDepth]) ;
        
        % remove abdominal segment?
        if removeASFlag
           imMat(AS_row_start:end, :, :) = 0 ;  
        end
        % call a function to do the binarization
        bwMat = binarize_confocal_stack(imMat) ;
    else
        % -------------------------------------------
        %% otherwise we don't need to worry as much about binarization
        % ----------------------------------------
        % convert to 8-bit
%         if isa(imMat,'uint8')
%             pixelScaleFactor = 1 ;
%         elseif isa(imMat,'uint16')
%             pixelScaleFactor = 1/256 ; 
%         end
%         imMat = uint8(pixScaleFactor.*double(imMat)) ; 
        imMat = uint8(255*mat2gray(imMat)) ; 
        
        % -----------------------------------------------------------------
        % resize to new dimensions (including accounting for non-cubic 
        %  voxels)
        imMat = imresize3(imMat,[imHeight, imWidth, imDepth]) ;
        
        % -----------------------------------------------------------------
        % do some extra processing?
        if extraProcessFlag
            imMat = imtophat(imMat, strel('sphere',4)) ;
            imMat = imadjustn(imMat,[0.0, 0.5]) ;
%             test = imsubtract(imadd(imMat,imtophat(imMat,strel('sphere',4))),...
%                 imbothat(imMat,strel('sphere',4)));
            %maxLevelCurr = 0.0675 ;
        else
            %maxLevelCurr = maxLevel ;
        end
        % -----------------------------------------------------------------
        %% convert to binary
        otsuLevel = graythresh(imMat) ;
        level = min([otsuLevel, maxLevel]) ;
        bwMat = imbinarize(imMat,level) ;
        
        if extraProcessFlag
           bwMat = imclose(bwMat, strel('sphere',2)) ; 
           bwMat = bwmorph3(bwMat,'fill') ;
        end
        
        % clean image a little
        bwMat = bwareaopen(bwMat, minObjSize) ;
        
%         % make sure we didn't leave any gaps
%         bwMat = bwmorph3(bwMat,'fill') ; 
        
        % make sure we're not removing too many pixels
        if (sum(bwMat(:)) < 1)
            % if we've opened the image and have no pixels left, just take
            % largest object
           bwMat = imbinarize(imMat,level) ;
           CC = bwconncomp(bwMat) ; 
           CCSizes = arrayfun(@(x) length(x.PixelIdxList), CC) ;
           maxCCSize = max(CCSizes) ; 
           
           % open image to take just the largest connected component (so
           % we're still cleaning up any fragments)
           bwMat = bwareaopen(bwMat, maxCCSize - 1) ; 
        end
        
    end
    
    %% filter out objects that don't connect to neck?
    if DNFlag
       CC = bwconncomp(bwMat) ; 
       bw_props = regionprops3(CC, 'BoundingBox','VoxelIdxList') ; 
       
       boundingBoxHeights = bw_props.BoundingBox(:,2) ; 
       good_idx = (boundingBoxHeights < DN_top_dist_thresh) ; 
       
       bwMat = false(size(bwMat)) ; 
       bwMat(cell2mat(bw_props.VoxelIdxList(good_idx))) = true ; 
       
    end
    
    
    % -------------------------------------------------------------------
    %% make symmetrical about long body axis?
    if symFlag
        % check if image already seems bilateral
        symCheck = sum(bwMat & flip(bwMat,2),'all')/sum(bwMat(:)) ;
        %disp(symCheck)
        
        % if it seems unilateral, symmetrize
        if symCheck < symThresh
            bwMat = bwMat | flip(bwMat,2) ;
        end
        
    end
    
    % -----------------------------------------------------
    %% save results (and plot?)
    if debugFlag
        % initialize figure
        figUnits = 'inches' ; 
        figPosition = [5.4583, 1.8333, 6.2917, 7.8125] ;
        h_IN = figure('PaperPositionMode','auto','MenuBar','none',...
            'ToolBar','none','DockControls','off','Units',figUnits,...
            'OuterPosition',figPosition) ;

        % initialize axis
        ax = gca ;
        hold(ax,'on')

        % create a parent hgtform 
        parent = hgtransform('Parent',ax);
        
        % plot params
        vncAlpha = 0.05 ;
        neuronAlpha = 1.0 ; 
        inputColor = lines(1) ;
        imSize = size(bwMat) ; 
        
        % draw VNC background and neuron
        drawVNCBackgroundToAxis(ax, h_IN, [], [], vncAlpha, parent) ;
        ax = drawNeuronToAxis(ax, bwMat, inputColor, neuronAlpha, ...
            imSize, parent) ;
        
        exportgraphics(h_IN, fullfile(bwPath, [fn '_bw_plot.png']),...
                'Resolution', 500)
        close all
    end
    save(bwSavePath,'bwMat')
   
    toc
end





