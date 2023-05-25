% -------------------------------------------------------------------------
% script to generate plot binarized neuron over a vnc template to give a
% feel for where things are

% -------------------------------------------------------------------------
%% data path info
%dataRoot = 'D:\Fly Imaging\Erica Interneuron Stacks\IN\' ;
dataRoot = 'D:\Fly Imaging\Erica Interneuron Stacks\DN\' ;
savePath = 'D:\Fly Imaging\Erica Interneuron Stacks\DN\binary_figs\' ; 
stockListPath = 'D:\Fly Imaging\Erica Interneuron Stacks\DN\' ; 
stockListFilename = 'candidate_dn_lines.csv' ;

tiffDir = dir([dataRoot '*.tif']) ;
nrrdDir = dir([dataRoot '*.nrrd']) ;
imDir = vertcat(nrrdDir, tiffDir) ;
N_im = length(imDir) ;

% save over? plot?
overWriteFlag = true ;

% image info
N_channels = 3 ; % the segmented images are just one channel (green)
nc82_channel_ind = 1 ;
    
% processing params
voxDimDefault = [0.461, 0.461, 0.7] ;
shrinkFactor = 7 ; % how much to reduce volume by 

% -------------------
% plot params
az_dorsal = -180 ; 
el_dorsal = 90 ;

az_coronal = 0 ; 
el_coronal = 0 ; 

vncColor = 0.7*[1, 1, 1] ; 
vncAlpha = 0.2 ; 

neuronColor = [0, 0, 0] ; 
neuronAlpha = 1.0 ; 

figPosition = [1921, 41, 1920, 963] ; 

% directory for storing binarize/resized images
bwPath = fullfile(dataRoot, 'binarized_new') ;
if ~exist(bwPath, 'dir')
    keyboard
end
bwDir = dir(fullfile(bwPath, '*_bw.mat')) ; 

% ----------------------------------------------------------
%% load catalog to get list of good images (or just do all)
stockListPathFull = fullfile(stockListPath, stockListFilename ) ; 

if exist(stockListPathFull,'file')
    % import stock list data
    opts = detectImportOptions(stockListPathFull);
    stockListData = readtable(stockListPathFull, opts) ;
    
    % get names of drivers for which we have decent images
    driverNames = stockListData.LineName ; 
    goodImageIdx = logical(stockListData.goodImage) ; 
    goodDrivers = driverNames(goodImageIdx) ; 
    
    % find corresponding indices in directory structure
    bw_idx = false(1, length(bwDir)) ; 
    for j = 1:length(bwDir)
        % get driver name from filename (for bw mat)
        fn_curr = bwDir(j).name ;
        [ind1, ind2] = regexp(fn_curr, 'SS\d\d\d\d\d') ; 
        if isempty(ind1)
            continue
        end
        driver_curr = fn_curr(ind1:ind2) ;
        
        % see if driver is in list of "good drivers"
        match_idx = cellfun(@(y) strcmp(y, driver_curr), goodDrivers) ; 
        
        if (sum(match_idx) > 0)
           bw_idx(j) = true ;  
        end
    end
else
    bw_idx = true(1,length(bwDir)) ; 
end
% -----------------------------------------------------------
%% loop over images and process them
bw_ind = find(bw_idx) ; 

for i = 6  %1:length(bw_ind) %N_im %1:N_im
    tic
    ind = bw_ind(i) ;
    fprintf('Processing %d / %d stacks \n', i, length(bw_ind))
    
    % ----------------------------------
    %% load data
    % get filename for image
    dataFilename = fullfile(imDir(ind).folder, imDir(ind).name) ;
    [~, fn, ext] = fileparts(dataFilename) ;
    % get filename for bw mat
    bwDataPath = fullfile(bwPath, [fn '_bw.mat']) ;
    % get filename to save to 
    savePathDorsal = fullfile(savePath, [fn '_seg_attempt_dorsal.png']) ; 
    savePathCoronal = fullfile(savePath, [fn '_seg_attempt_coronal.png']) ; 
    
    % also get name of current driver and DN name
    [ind1, ind2] = regexp(fn, 'SS\d\d\d\d\d') ; 
    driver = fn(ind1:ind2) ; 
    
    fn_split = strsplit(fn(1:(ind1-2)),'_') ; 
    line_name = strjoin(fn_split, ', ') ; 
    
    if (exist(savePathDorsal,'file') || exist(savePathCoronal,'file')) && ...
            ~overWriteFlag
       fprintf('Already created figure for %s \n', bwDataPath)
       continue
    end
    
    % load binarized neuron 
    neuronBW = importdata(bwDataPath) ; 
    
    % get image dimensions (just in case)
    [imHeight, imWidth, imDepth] = size(neuronBW) ; 
    % ---------------------------------------------------------------
    % if it's a multi-channel image, grab vnc channel. otherwise use
    % pre-saved template
    if N_channels > 1
        % load image
        % load image
        if strcmp(ext,'.tif')
            imMat  = myReadTiff(dataFilename, N_channels) ;
        elseif strcmp(ext,'.nrrd')
            [imMat, ~] = nrrdread(dataFilename) ;
        else
            fprintf('Invalid file extension: %s \n', ext)
            continue
        end
        
        % get only nc82 channel
        imMat = squeeze(imMat(:,:,nc82_channel_ind,:)) ;
        
        % resize image
        imMat = imresize3(imMat,[imHeight, imWidth, imDepth]) ;
        
        % binarize 
        vncBW = binarize_confocal_stack(imMat) ; 
        
    else
        fprintf('Under construction! \n')
        keyboard
    end
    
    % ---------------------------------------------------------------
    %% plot vnc and neuron
    h_main = figure('Position',figPosition) ;
    hold on
    ax = gca ; 
    ax = drawNeuronsAndVNC(ax, h_main, vncBW, {neuronBW}, ...
        neuronColor, vncColor, vncAlpha) ;
    
    % add driver label
    xlim = get(ax, 'xlim') ; 
    ylim = get(ax, 'ylim') ;
    zlim = get(ax, 'zlim') ; 
    message = sprintf([line_name '\n' driver]) ;
    text(xlim(2), ylim(1), zlim(2), message,'fontSize',16) 
    % -----------------------------------------------------------------
    %% save results
    print(h_main,savePathDorsal,'-dpng','-r500') 
    view(az_coronal, el_coronal)
    print(h_main,savePathCoronal,'-dpng','-r500') 
    
    toc
    fprintf('Completed %d / %d stacks \n', i, length(bw_ind))
    close all
end