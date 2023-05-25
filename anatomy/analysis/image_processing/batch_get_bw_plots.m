% -------------------------------------------------------------------------
% script to generate new png plots for binarized images in a given folder
% -------------------------------------------------------------------------
%% data types and params
dataTypeList = {'IN'} ;  %{'MN', 'IN', 'DN', 'VUM'} ; % 'IN'

symFlag = true ;
% save over? plot?
overWriteFlag = false ;

% plot params:
figUnits = 'inches' ;
figPosition = [5.4583, 1.8333, 6.2917, 7.8125] ;
vncAlpha = 0.05 ;
neuronAlpha = 1.0 ;
inputColor = lines(1) ;

% initialize figure window
h_IN = figure('PaperPositionMode','auto','MenuBar','none',...
    'ToolBar','none','DockControls','off','Units',figUnits,...
    'OuterPosition',figPosition) ;

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
    %% loop over images and make plots
    for ind = 1:N_bw
        tic
        
        fprintf('Processing %d / %d stacks \n', ind, N_bw)
        % get filename for image
        dataFilename = fullfile(bwDir(ind).folder, bwDir(ind).name) ;
        [~, fn, ext] = fileparts(dataFilename) ;
        % set save path
        plotSavePath = fullfile(savePath, [fn '_bw_plot.png']) ;
        
        % check if we've already converted this one
        if (~overWriteFlag) && exist(plotSavePath,'file')
            fprintf('Already have plot for: %s \n', dataFilename)
            continue
        end
        % -----------------------------------------------------------
        % load image
        bwMat = importdata(dataFilename) ;
        
        % -------------------------------------------------------------------
        %% make plot 
        % initialize axis
        set(0, 'CurrentFigure', h_IN)
        ax = gca ;
        hold(ax,'on')

        % create a parent hgtform 
        parent = hgtransform('Parent',ax);
        
        % get size of current image
        imSize = size(bwMat) ; 
        
        % draw VNC background and neuron
        drawVNCBackgroundToAxis(ax, h_IN, [], [], vncAlpha, parent) ;
        ax = drawNeuronToAxis(ax, bwMat, inputColor, neuronAlpha, ...
            imSize, parent) ;

        % -------------------------------------------
        %% save plot and clear figure window
        exportgraphics(h_IN, plotSavePath, 'Resolution', 500)
        clf(h_IN)
        
        fprintf('Saved %d / %d figures \n', ind, N_bw)
        toc
    end
end
