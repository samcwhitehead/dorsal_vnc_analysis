% -------------------------------------------------------------------------
% script to generate mips of images in a directory -- writing this
% separately so we can generate the images and have them stored for later
% use (rather than loading a stack every time one is needed, which is slow)
% -------------------------------------------------------------------------
dataRoot = 'D:\Fly Imaging\Erica Interneuron Stacks\IN\' ;
savePath = fullfile(dataRoot, 'MIP') ; 
if ~exist(savePath,'dir')
    mkdir(savePath) ;
end
%dataRoot = 'D:\Fly Imaging\Erica Interneuron Stacks\DN\' ;
tiffDir = dir([dataRoot '*.tif']) ;
nrrdDir = dir([dataRoot '*.nrrd']) ;
imDir = vertcat(nrrdDir, tiffDir) ;
N_im = length(imDir) ;

% what type of image to save to 
imFileExt = '.png' ; 

% save over? 
overWriteFlag = false ; 

% channels per image?
N_channels = 1 ; 

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
    mipSavePath = fullfile(savePath, [fn '_mip']) ;
    
    % check if we've already converted this one
    if (~overWriteFlag) && exist([mipSavePath, imFileExt],'file')
        fprintf('Already analyzed file: %s \n', dataFilename)
        continue       
    end
    % -----------------------------------------------------------
    % load image
    if strcmp(ext,'.tif')
        imMat  = myReadTiff(dataFilename, N_channels) ;
    elseif strcmp(ext,'.nrrd')
        [imMat, ~] = nrrdread(dataFilename) ;
    else
        fprintf('Invalid file extension: %s \n', ext)
        continue
    end
    
    % ---------------------------------------------------------------
    % process image
    imMat = uint8(255*mat2gray(imMat)) ; % convert to 8-bit
    imMat = imresize3(imMat,[imHeight, imWidth, imDepth]) ;  % resize
    
    % get mip 
    imMIP = max(imMat,[],3) ; 
    
    % ---------------------------------------------------------------
    % save mip
    imwrite(imMIP, [mipSavePath, imFileExt])
    save([mipSavePath, '_data.mat'], 'imMIP')
    
    toc
end
