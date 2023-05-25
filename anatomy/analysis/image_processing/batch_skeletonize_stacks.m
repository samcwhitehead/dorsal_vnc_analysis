% -------------------------------------------------------------------------
% script to loop through segmented interneuron stacks from Erica, binarize
% them, resize them, and then save the logical arrays. the hope is that
% this will speed up overlap computations, since much of the overhead seems
% to be i/o
% -------------------------------------------------------------------------
%% data path info
dataType = 'MN' ; % 'IN'
dataRoot = ['D:\Fly Imaging\Erica Interneuron Stacks\' dataType] ;

%dataRoot = 'D:\Fly Imaging\Erica Interneuron Stacks\DN\' ;
symFlag = true ; 

if symFlag
    dataPath = fullfile(dataRoot, 'binarized_new') ; 
    savePath = fullfile(dataRoot, 'skeletonized_sym') ; 
else
    dataPath = fullfile(dataRoot, 'binarized_non_sym') ; 
    savePath = fullfile(dataRoot, 'skeletonized_non_sym') ; 
end
if ~exist(savePath,'dir')
    mkdir(savePath)
end
bwDir = dir(fullfile(dataPath, '*bw.mat')) ;
N_bw = length(bwDir) ;

% save over? plot?
overWriteFlag = true ;
debugFlag = false ;
plotFlag = true ;

if strcmp(dataType, 'DN')
    THR = 25 ; % minimum branch length; used to filter artifacts % WAS 10 (previously 0)
else
    THR = 10 ;  % for INs
end 
% ----------------------------------------------
% loop over images and skeletonize them
for ind = 1:N_bw
    tic
    
    fprintf('Processing %d / %d stacks \n', ind, N_bw)
    % get filename for image
    dataFilename = fullfile(bwDir(ind).folder, bwDir(ind).name) ;
    [~, fn, ext] = fileparts(dataFilename) ;
    % set save path
    skeletonSavePath = fullfile(savePath, [fn '_skel_struct.mat']) ;
    
    % check if we've already converted this one
    if (~overWriteFlag) && exist(skeletonSavePath,'file')
        fprintf('Already analyzed file: %s \n', dataFilename)
        continue
    end
    % -----------------------------------------------------------
    % load image
    bwMat = importdata(dataFilename) ;
    
    % -------------------------------------------------------------------
    % do skeletonization
    skel_struct = mySkelAndGraph(bwMat, THR, debugFlag) ;
    
    if debugFlag
        keyboard
    end
    
    % -------------------------------------------
    % save output
    save(skeletonSavePath,'skel_struct')
    
    % -------------------------------------------
    % make plot to check skeletonization process?
    if plotFlag
       h_temp = figure ; 
       drawSkelGraph([],skel_struct) ;
       title(fn,'Interpreter','none') 
       saveFigPath = fullfile(savePath, [fn '_skel_plot.png']) ;
       print(h_temp, saveFigPath, '-dpng', '-r300') 
       close(h_temp)
    end
    toc
end






