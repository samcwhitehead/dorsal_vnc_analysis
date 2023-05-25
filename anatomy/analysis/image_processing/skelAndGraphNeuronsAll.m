% -------------------------------------------------------------------------
% quick and dirty script to loop through binarized neuron images and get
% skeleton/graph structure
% -------------------------------------------------------------------------
rootPath = 'D:\Fly Imaging\Erica Interneuron Stacks\' ; 
neuronType = 'IN' ; 
dataRoot = fullfile(rootPath, neuronType) ; 
bwDir = dir(fullfile(dataRoot,'binarized_new', '*bw.mat')) ;

saveFlag = true ;
plotFlag = true ; 
overWriteFlag = false ; 

THR = 5 ; % minimum skeleton branch length
N_stacks = length(bwDir) ; 

% --------------------------------------
% loop through bw images and skeletonize
for ind = 1:N_stacks
    tic
    fprintf('Processing %d / %d stacks \n', ind, N_stacks)
    
    % -------------------------------------------------------------
    % get filename for image
    bwPath = bwDir(ind).folder ; 
    dataFilename = fullfile(bwPath, bwDir(ind).name) ;
    [~, fn, ~] = fileparts(dataFilename) ; 
    
    fn_split = strsplit(fn,'_') ; 
    fn_new = strjoin(fn_split(1:end-1),'_') ;
    % set save path
    bwSavePath = fullfile(bwPath, [fn_new '_skel_struct.mat']) ; 
    if exist(bwSavePath,'file') && ~overWriteFlag
       fprintf('Already analyzed %s -- skipping \n', fn_new) 
    end
    % -------------------------------------------------------------
    % load data
    bwMat = importdata(dataFilename) ; 

    % -------------------------------------------------------------
    % compute skeleton and graph representation
    [skel_struct, h_main] = mySkelAndGraph(bwMat, THR, plotFlag) ; 
    
    % -------------------------------------------------------------
    % save results?
    % figure:
    if plotFlag 
        if saveFlag 
            figSavePath = fullfile(bwPath, [fn_new '_skel_plot.fig']) ; 
            savefig(h_main, figSavePath)
        end
        close(h_main)
    end
    
    % data:
    if saveFlag
       save(bwSavePath, 'skel_struct')
    end
    toc
end