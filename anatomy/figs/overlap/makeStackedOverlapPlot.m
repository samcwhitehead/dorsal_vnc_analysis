% -------------------------------------------------------------------------
% script to generate stackable overlap plots using combined overlap tables
% (e.g. [MN and IN] vs [DN])
%
% TO DO:
%   - get neuron-type-specific sorting indices from full list
%   - figure out how to set consistent data aspect ratio
%
% -------------------------------------------------------------------------
% path info
% path info
[mfilePath, ~, ~] = fileparts(mfilename('fullpath')) ; 
figDirectory = fileparts(mfilePath) ;
parentDirectory = fileparts(figDirectory) ; 
dataRoot = fullfile(parentDirectory, 'data') ;
overlapPath = fullfile(dataRoot, 'overlap_calculations') ;

savePath = fullfile(mfilePath, 'output', 'stacked_overlap') ;
if ~exist(savePath, 'dir')
    mkdir(savePath)
end

% which combinations of overlap matrices to plot
dataNameList1 = {'IN', 'MN', 'VUM'} ; 
dataNameList2 = {'IN', 'MN', 'VUM'} ; 
% dataNameList1 = {'IN', 'MN', 'VUM'} ; 
% dataNameList2 = {'DN'} ; 

% params
saveFlag = true ; % save plots
normType = 'none' ; % how to normalize overlap data? % 'sym_size' ; | 'none'
maskedFlag = true ; % use masked overlap?
threshVal = 300 ; % minimum value of overlap. 500 ~ a 3.7x3.7x3.7 micron cube
globalSortFlag = true ; % sort rows/cols consistently across blocks?
rmColorBarFlag = true ; 

% some plot params
colorTickLabelsFlag = false ; 
clim = [0, 1.65e4] ; % [] ; 
axisFontSizeSmall = 4.5 ; 
if any(strcmp(dataNameList1, dataNameList2))
    data_aspect_ratio = [2,2,1] ; 
else
    data_aspect_ratio = [2,3,1] ; % will maybe depend on neuron types?
end
% search expression for overlap filenames
searchStr = '(?<dataName1>[A-Z]{2,3})_(?<dataName2>[A-Z]{2,3})_overlap_table' ;
% update overlap path if we're using masked overlap (as we probably should
% be)
if maskedFlag
   overlapPath = fullfile(overlapPath, 'masked') ;  
end

% filename to save to (don't use save option within function)
saveFnStr = '%s_%s_%s_overlap' ;

% base figure position and inches of plot per neuron
figPosition = [0.5, 0.5, 6.5, 11.0] ;
% inchPerNeuronX = 0.1 ; 
% inchPerNeuronY = 0.1 ; 

% make save directory if it doesn't exist already. 
% NB: we want to make a folder for each group of plots that will go into
% one set of stacked plots
saveFolder = sprintf('%s_%s_overlap',strjoin(dataNameList1,'-'), ...
    strjoin(dataNameList2,'-')) ;  %
savePath = fullfile(savePath, saveFolder) ; 
if ~exist(savePath,'dir')
    mkdir(savePath)
end

% -----------------------------------------------------------------------
%% load neuron color struct, if using
if colorTickLabelsFlag
    % load neuron color struct struct corresponding to neuron TYPE
    colorStructPath = fullfile(dataRoot, 'neuronTypeColorStruct.mat') ;
    neuronColorStruct = importdata(colorStructPath) ;
end

% -----------------------------------------------------------------------
%% get number of neurons per type
% will be used for applying sorting indices
N_neurons1 = zeros(length(dataNameList1)+1,1) ; % adding one to make looping easier
for q = 1:length(dataNameList1)
    dataPathCurr = fullfile(dataRoot, 'vox_size_tables', dataNameList1{q}) ;
    % load vox size table
    vox_sz_table_curr = importdata(fullfile(dataPathCurr, ...
        'vox_size_table_masked.mat')) ; 
    % get length to determine neuron count
    N_neurons1(q+1) = height(vox_sz_table_curr) ; 
end

N_neurons2 = zeros(length(dataNameList2)+1,1) ; % adding one to make looping easier
for q = 1:length(dataNameList2)
    dataPathCurr = fullfile(dataRoot, 'vox_size_tables', dataNameList2{q}) ;
    % load vox size table
    vox_sz_table_curr = importdata(fullfile(dataPathCurr, ...
        'vox_size_table_masked.mat')) ; 
    % get length to determine neuron count
    N_neurons2(q+1) = height(vox_sz_table_curr) ; 
end

% use this to get x and y limits for plots
xlim = [0.5, max(N_neurons2) + 0.5] ; 
ylim = [0.5, max(N_neurons1) + 0.5] ; 

% -----------------------------------------------------------------------
%% get global sorting pattern (if using)
% okay, this is more complicated than i thought, esp. when we're doing both
% row and column blocks 
if globalSortFlag
    % call plot function just for sorting indices
    [hfig, ax, row_ind, col_ind] = makeNormOverlapPlotFunc(dataNameList1, ...
        dataNameList2, savePath, dataRoot, false, normType, maskedFlag, ...
        threshVal, false, clim, figPosition, 'new_cluster') ;
    
    % get color bar limits for combined plot
    if isempty(clim)
        clim = ax.CLim ; 
    end
    
    % close combined plot
    close(hfig)
    
    % set value of matSortCase for subsequent plots (which we'll save)
    matSortCase = 'manual' ; 
    
else
    % otherwise, just need to set sorting method. this should be
    % 'new_cluster' since we want each block to be clustered as well as
    % possible
    matSortCase = 'new_cluster' ;
    row_ind = [] ; 
    col_ind = [] ; 
end

% -----------------------------------------------------------------------
%% loop over neuron combinations to make figure
% NB: as in "loadCombinedOverlapMat.m" we're going to work with stack
% columns as opposed to looping through the full grid at once

% get possible combinations of neuron types
[xx, yy] = meshgrid(1:length(dataNameList1), 1:length(dataNameList2)) ;
comb_ind = [xx(:), yy(:)] ;

% number of combinations
N_comb = size(comb_ind,1) ;

% number of column blocks
N_col_blocks = length(unique(comb_ind(:,2))) ;

% set width for the figures within each column block based on max # neurons
figPosition(3) = 0.1226*max(N_neurons2) ;
    
% start by looping over all combinations that involve the same dataName2
for m = 1:N_col_blocks
    % find indices in comb_ind that correspond to current dataName2
    ind_curr = find(comb_ind(:,2) == m) ; 
    
    % get dataName that corresponds to this column
    dataName2 = dataNameList2{m} ; 
    
    % initialize storage for graphics handles of each matrix within column
    % block 
    h_array = gobjects(length(ind_curr),1) ; 
    ax_array = gobjects(length(ind_curr),1) ; 
    
    
    
    % get current sorting indices for column
    if globalSortFlag
        % first get indices of column neurons in unsorted matrix
        col_neuron_nums = (1:N_neurons2(m+1)) + sum(N_neurons2(1:m)) ; 
        
        % then pull out the sorting ind for these neurons from the
        % global sort list
        neuron_col_idx = (col_ind >= col_neuron_nums(1)) & ...
            (col_ind <= col_neuron_nums(end)) ;
        col_ind_curr = col_ind(neuron_col_idx) ;
        
        % subtract off constant value so all of our col indices go from
        % 1 to N (where N is number of current column neurons)
        col_ind_curr = col_ind_curr - min(col_ind_curr) + 1 ;
        
    else
       % otherwise, leave empty, since w/o global sort each matrix gets its 
       % own clustering (calculated within plot function)
       col_ind_curr = [] ; 
    end
    
    % -------------------------------------------------
    % loop over rows within column block
    for n = 1:length(ind_curr)
        % read out current row data name
        dataName1 = dataNameList1{comb_ind(ind_curr(n),1)} ; 
        
        % get appropriate sorting index, if using global sort
        if globalSortFlag
            % first get indices for current row neurons in unsorted matrix 
            row_neuron_nums = (1:N_neurons1(n+1)) + sum(N_neurons1(1:n)) ; 
            
            % then pull out the sorting ind for these neurons from the
            % global sort list
            neuron_row_idx = (row_ind >= row_neuron_nums(1)) & ...
                (row_ind <= row_neuron_nums(end)) ;
            row_ind_curr = row_ind(neuron_row_idx) ; 
            
            % subtract off constant value so all of our row indices go from
            % 1 to N (where N is number of current row neurons)
            row_ind_curr = row_ind_curr - min(row_ind_curr) + 1 ; 
            
        else
            % otherwise we can keep row_ind_curr empty, since, if we're not
            % using global sorting, we'll let the matrix cluster itself
            row_ind_curr = [] ; 
        end
        
        % ----------------------------------------------------
        %% make plot
        [h_array(n), ax_array(n)] = makeNormOverlapPlotFunc(dataName1, ...
            dataName2, savePath, dataRoot, false, normType, maskedFlag, ...
            threshVal, false, clim, figPosition, matSortCase, ...
            row_ind_curr, col_ind_curr) ;
        
        % turn off "axis tight"
        axis(ax_array(n),'normal')
        
        % try to color tick labels by neuron type? NB: different from method in
        % function, which tries to use cluster coloring
        if colorTickLabelsFlag
            % color tick labels
            ax_array(n) = colorNeuronTickLabels(ax_array(n), ...
                neuronColorStruct) ;
        end
        
        % try to increase height of axis to fit names better?
%         if all(strcmp(dataName1,dataName2))
        ax_pos = get(ax_array(n),'Position') ;
        ax_pos(4) = ax_pos(4) + 0.065 ;
        set(ax_array(n),'Position',ax_pos)
        
        % also make square if we need to (same data types on X and Y)
        if any(strcmp(dataNameList1,dataNameList2))
            axis(ax_array(n),'equal') 
        end
        
        % get current x/y limits before changing them
        xlim_curr = get(ax_array(n),'xlim') ;
        ylim_curr = get(ax_array(n),'ylim') ; 
        
        % set axis limits
        set(ax_array(n), 'xlim', xlim, 'ylim', ylim)
        
        % adjust font size (this happens a little in function, but since we
        % may have an effectively larger number of neurons than is given in
        % each function call, need to adjust again)
        ax_array(n).XAxis.FontSize = axisFontSizeSmall ;
        ax_array(n).YAxis.FontSize = axisFontSizeSmall ;
        
        % move title so that it more closely abuts plot
        h_tit = get(ax_array(n),'Title') ; 
        tit_str = get(h_tit, 'String') ; 
        tit_str_split = strsplit(tit_str, ' ') ; 
        tit_str = strjoin(tit_str_split(1:(end-1)), ' ') ; 
        set(h_tit,'String', tit_str) ;
        tit_pos = get(h_tit, 'Position') ; 
        tit_pos(2) = ylim_curr(2) + 0.4 ; 
        tit_pos(1) = mean(xlim_curr) ; 
        set(h_tit, 'Position', tit_pos) ; 
        
        % remove colorbar?
        if rmColorBarFlag
            c = findobj(h_array(n), 'Type', 'colorbar') ; 
            set(c,'visible','off')
        end
        
    end
    
    % ------------------------------------------------------------------
    %% save column block plots
    % now want to save column block. to do this, we'll first need to set
    % data aspect ratio to be consistent across plots
%     data_aspect_list = zeros(length(ax_array), 3) ;
%     ax_pos_list = zeros(length(ax_array), 4) ;
%     equal_ax_pos = [] ; % initialize variable for position of square overlap matrix
%     % pb_aspect_list = zeros(length(ax_array), 3) ; 
%     for p = 1:length(ax_array)
%         data_aspect_list(p,:) = daspect(ax_array(p)) ; 
%         ax_pos_list(p,:) = get(ax_array(p),'Position') ; 
%         
%         dataName1 = dataNameList1{comb_ind(ind_curr(p),1)} ; 
%         if strcmp(dataName1,dataName2)
%            equal_ax_pos = ax_pos_list(p,:) ;  
%         end
%         % pb_aspect_list(p,:) = pbaspect(ax_array(p)) ; 
%     end
% %     
% %     % first sort out axis width
% %     if ~isempty(equal_ax_pos)
% %        for p = 1:length(ax_array)
% %           ax_pos_new = ax_pos_list(p,:) ; 
% %           ax_pos_new(3) = equal_ax_pos(3) ; 
% %           set(ax_array(p),'Position', ax_pos_new) ; 
% %        end
% %     end
% %     % get maximum data aspect ratio
% %     [~, max_ind] = max(data_aspect_list(:,2)./data_aspect_list(:,1)) ; 
% %     

    % set all axes to have same data aspect ratio
    %NB: for some reason, this messes up axis equal plots
    for p = 1:length(ax_array)
        %         daspect(ax_array(p), data_aspect_list(max_ind,:)) ;
        daspect(ax_array(p), data_aspect_ratio) ;
    end
    
    % save results
    if saveFlag
        % suffix string (tells whether or not we applied global sorting)
         if globalSortFlag
            suffixStr = '_globalSort' ;
         else
             suffixStr = '' ; 
         end
         
         % loop over figure windows
         for q = 1:length(h_array)
             % get current save filename
             dataName1 = dataNameList1{comb_ind(ind_curr(q),1)} ; 
             saveFn = sprintf(saveFnStr, dataName1, dataName2, normType) ;
             
             % save FIG file
             savefig(h_array(q), fullfile(savePath, ...
                 [saveFn suffixStr '.fig'])) ;
             
             % save PNG file
             exportgraphics(h_array(q), fullfile(savePath, ...
                 [saveFn suffixStr '.png']),'Resolution',600) ;
             
             % save JPEG file
             exportgraphics(h_array(q), fullfile(savePath, ...
                 [saveFn suffixStr '.jpg']),'Resolution',600) ;
         end
        
         % close open figure windows
         close(h_array(:))
    end
    
end


