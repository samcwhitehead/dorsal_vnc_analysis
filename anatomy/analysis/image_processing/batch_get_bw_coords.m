% -------------------------------------------------------------------------
% script to loop through binarized interneuron stacks from Erica and get a
% list of coordinates of "on" voxels. this will reduce the data file size
% for each neuron, and hopefully speed up overlap calculation, since much
% of the overhead seems to be i/o
% -------------------------------------------------------------------------
%% options
dataTypeList = {'VUM'} ;  %{'MN', 'IN', 'DN', 'VUM'} ; % 'IN'

symFlag = true ;
% save over? plot?
overWriteFlag = false ;

% ------------------------------------------------------------------
%% loop over data types
for d = 1:length(dataTypeList)
    
    % get path info for current data type
    dataType = dataTypeList{d}  ;
    dataRoot = ['D:\Fly Imaging\Erica Interneuron Stacks\' dataType] ;
    
    if symFlag || strcmpi(dataType,'VUM')
        dataPath = fullfile(dataRoot, 'binarized_new') ;
    else
        dataPath = fullfile(dataRoot, 'binarized_non_sym') ;
    end
    savePath = dataPath ; 
    
    % get directory of bw images for current data type
    bwDir = dir(fullfile(dataPath, '*bw.mat')) ;
    N_bw = length(bwDir) ;
    
    
    % ----------------------------------------------
    % loop over images and get their coordinates
    for ind = 1:N_bw
        tic
        
        fprintf('Processing %d / %d stacks \n', ind, N_bw)
        % get filename for image
        dataFilename = fullfile(bwDir(ind).folder, bwDir(ind).name) ;
        [~, fn, ext] = fileparts(dataFilename) ;
        % set save path
        coordsSavePath = fullfile(savePath, [fn '_coords.mat']) ;
        
        % check if we've already converted this one
        if (~overWriteFlag) && exist(coordsSavePath,'file')
            fprintf('Already analyzed file: %s \n', dataFilename)
            continue
        end
        % -----------------------------------------------------------
        % load image
        bwMat = importdata(dataFilename) ;
        
        % -------------------------------------------------------------------
        % get coordinates
        coords = find(bwMat(:)) ;
        
        % -------------------------------------------
        % save output
        save(coordsSavePath,'coords')
        
        toc
    end
end





