% -------------------------------------------------------------------------
% script to generate overlap plots using combined overlap tables (e.g. [MN
% and IN] vs [DN])
% -------------------------------------------------------------------------
% path info
dataRoot = 'D:\Fly Imaging\Erica Interneuron Stacks\' ;
overlapPath = fullfile(dataRoot, 'overlap_calculations') ;
savePath = ['D:\Dropbox\Paper Manuscripts\Janelia vnc paper\'...
    'figure drafts\DN_overlap_fig\combined_overlap_matrices\'] ;

% which combinations of overlap matrices to plot
% dataNames1 = {'IN', 'MN', 'VUM'} ; 
% dataNames2 = {'IN', 'MN', 'VUM'} ; 
dataNames1 = {'IN', 'MN', 'VUM'} ; 
dataNames2 = {'DN'} ; 

% which plot case to use. we either:
%   1) plot combined overlap data as one big matrix
%   2) make stacked overlap plots with a common column sorting
%   3) make stacked overlap plots each with their own sorting

% TEMP test the way we sort rows/columns
testSortFlag = false ; 
if testSortFlag
   savePath = fullfile(savePath, 'sort_test') ;  
end

% params
saveFlag = false ; % save plots
normType = 'none' ; % how to normalize overlap data? % 'sym_size' ; | 'none'
maskedFlag = true ; % use masked overlap?
threshVal = 300 ; % minimum value of overlap. 500 ~ a 3.7x3.7x3.7 micron cube
matSortCase = 'new_cluster' ; % how to sort matrix rows/cols

% some plot params
colorTickLabelsFlag = false ; 
clim = [] ; 

% search expression for overlap filenames
searchStr = '(?<dataName1>[A-Z]{2,3})_(?<dataName2>[A-Z]{2,3})_overlap_table' ;
% update overlap path if we're using masked overlap (as we probably should
% be)
if maskedFlag
   overlapPath = fullfile(overlapPath, 'masked') ;  
end

% filename to save to (don't use save option within function)
saveFn = sprintf('%s_%s_%s_overlap',strjoin(dataNames1,'-'), ...
    strjoin(dataNames2,'-'), normType) ;  %

% base figure position and inches of plot per neuron
figPosition = [0.5, 0.5, 6.5, 11.0] ;
inchPerNeuronX = 0.1 ; 
inchPerNeuronY = 0.1 ; 

% make save directory if it doesn't exist already
if ~exist(savePath,'dir')
    mkdir(savePath)
end
% -----------------------------------------------------------------------
%% make figure
[h_overlap, ax] = makeNormOverlapPlotFunc(dataNames1, dataNames2, ...
       savePath, dataRoot, false, normType, maskedFlag, threshVal, ...
       false, clim, figPosition, matSortCase) ;
   
% try to color tick labels by neuron type? NB: different from method in
% function, which tries to use cluster coloring
if colorTickLabelsFlag
    % load neuron color struct struct corresponding to neuron TYPE
    colorStructPath = fullfile(dataRoot, 'neuronTypeColorStruct.mat') ; 
    neuronColorStruct = importdata(colorStructPath) ; 
    
    % color tick labels
    ax = colorNeuronTickLabels(ax, neuronColorStruct) ;
end

% ------------------------------------------------------------------------
%% save plot?
% combined figure will probably need to be bigger than screen size, so make
% figure invisible then resize and save
if saveFlag
    %  make figure invisible
    set(h_overlap, 'visible', 'off') ; 
    
    % try to increase height of axis
    if all(strcmp(dataNames1,dataNames2)) 
        axis(ax,'equal')
    else
        ax_pos = get(ax,'Position') ;
        ax_pos(4) = ax_pos(4) + 0.065 ;
        set(ax,'Position',ax_pos)
    end
    
    % get new figure size
    n_rows = length(get(ax,'YTick')) ; 
    n_cols = length(get(ax,'XTick')) ; 
    
    figPositionNew = figPosition ; 
    figPositionNew(3) = inchPerNeuronX*n_cols ; 
    figPositionNew(4) = inchPerNeuronY*n_rows ; 
    
    set(h_overlap, 'OuterPosition', figPositionNew) ; %, ...
        %'PaperPosition', figPositionNew) ; 
    
    % if we're increasing figure size dramatically, try to reduce font
    if (figPositionNew(3) > 15)
       ax.XAxis.FontSize = 3.0 ;
    end
    if (figPositionNew(4) > 15)
       ax.YAxis.FontSize = 3.0 ;
    end
    
%     set(h_overlap, 'PaperPositionMode','manual',...
%         'PaperPosition',figPositionNew) ; 
    
%     % try to resize axis?
%     if all(strcmp(dataNames1, dataNames2))
%         axis equal
%     else
%         pause(1)
%         ax_pos = get(ax,'Position') ;
%         ax_pos(4) = ax_pos(4) + 0.065 ;
%         set(ax,'Position',ax_pos)
%     end
    
    % save output
    suffixStr = '_sort' ;
%     if ~testSortFlag
        savefig(h_overlap, fullfile(savePath, [saveFn suffixStr '.fig'])) ;
        exportgraphics(h_overlap, fullfile(savePath, ...
            [saveFn suffixStr '.png']),'Resolution',600) ;
%     end
    exportgraphics(h_overlap, fullfile(savePath, ...
        [saveFn suffixStr '.jpg']),'Resolution',600) ;
end