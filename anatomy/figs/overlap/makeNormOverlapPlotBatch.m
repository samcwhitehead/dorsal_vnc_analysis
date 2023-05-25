% -------------------------------------------------------------------------
% script to run "makeNormOverlapPlotFunc" on various combinations of neuron
% types
% -------------------------------------------------------------------------
% path info
[mfilePath, ~, ~] = fileparts(mfilename('fullpath')) ; 
figDirectory = fileparts(mfilePath) ;
parentDirectory = fileparts(figDirectory) ; 
dataRoot = fullfile(parentDirectory, 'data') ;
overlapPath = fullfile(dataRoot, 'overlap_calculations') ;

savePath = fullfile(mfilePath, 'output', 'overlap_plots') ;
if ~exist(savePath, 'dir')
    mkdir(savePath)
end

% params
saveFlag = true ; % save plots
normType = 'none' ; % how to normalize overlap data? % 'sym_size' ; | 'none'
maskedFlag = true ; % use masked overlap?
threshVal = 300 ; % minimum value of overlap. 500 ~ a 3.7x3.7x3.7 micron cube
matSortCase = 'new_cluster' ; % how to sort matrix rows/cols
figPositionAutoFlag = true ; % overrride other inputs and try to scale figure based on neuron count?
colorTickLabelsFlag = false ; 
clim = [0, 1.65e4] ; 

% if ~isempty(clim)
%    savePath = fullfile(savePath,'same_clim') ;  
% end

% search expression for overlap filenames
searchStr = '(?<dataName1>[A-Z]{2,3})_(?<dataName2>[A-Z]{2,3})_overlap_table' ;
% update overlap path if we're using masked overlap (as we probably should
% be)
if maskedFlag
   overlapPath = fullfile(overlapPath, 'masked') ;  
end

% base figure position 
figPosition = [0.5, 0.5, nan, 11.0] ; 

% ------------------------------------------------------------
% get info for the overlap calculations we have on file
overlapDir = dir(fullfile(overlapPath,'*overlap_table.mat')) ; 

% use regexp to find the two neuron types for each overlap file
fnNameMatches = arrayfun(@(x) regexp(x.name,searchStr,'names'), ...
    overlapDir, 'UniformOutput',0) ; 
dataNameList1 = cellfun(@(y) y.dataName1, fnNameMatches, ...
    'UniformOutput',0) ; 
dataNameList2 = cellfun(@(y) y.dataName2, fnNameMatches, ...
    'UniformOutput',0) ; 

if length(dataNameList1) ~= length(dataNameList2)
   fprintf('Error: different lengths for data name lists \n')
   keyboard
end

% --------------------------------------------------------------
% loop over data name combinations and make plots
for k = 1:length(dataNameList1)
   % read out current data types
   dataName1 = dataNameList1{k} ; 
   dataName2 = dataNameList2{k} ; 
   
   % get size for current figure
   % figure window size depends on size of neuron collection
   if figPositionAutoFlag
       figPosition = 'scaled' ;
   else
       
       switch dataName1  % FIGURE HEIGHT
           case 'DN'
               figPosition(4) = 6.5 ;
           case 'IN'
               figPosition(4) = 11.0 ;
           case 'MN'
               figPosition(4) = 4.5 ;
           case 'VUM'
               figPosition(4) = 5.5 ;
           otherwise
               fprintf('Invalid neuron type: %s \n', dataName1)
               keyboard
       end
       
       switch dataName2  % FIGURE WIDTH
           case 'DN'
               figPosition(3) = 6.5 ;
           case 'IN'
               figPosition(3) = 14.0 ;
           case 'MN'
               figPosition(3) = 3.25 ;
           case 'VUM'
               figPosition(3) = 5.5 ;
           otherwise
               fprintf('Invalid neuron type: %s \n', dataName2)
               keyboard
       end
   end
   
   % ------------------------------------------------------
   % make plot
   [h_overlap, ax] = makeNormOverlapPlotFunc(dataName1, dataName2, ...
       savePath, dataRoot, saveFlag, normType, maskedFlag, threshVal, ...
       colorTickLabelsFlag, clim, figPosition, matSortCase) ;
    
   % try to resize axis
   if all(strcmp(dataName1, dataName2))
       axis equal
   elseif ~strcmp(dataName1,'MN')
       pause(1)
       ax_pos = get(ax,'Position') ;
       ax_pos(4) = ax_pos(4) + 0.065 ;
       set(ax,'Position',ax_pos)
   end

   % close open figures
   if saveFlag
       close all
      %close(h_overlap) 
   end
   
   % print update
   fprintf('Completed plot for %s vs %s overlap (%d/%d)\n', dataName1, ...
       dataName2,k,length(dataNameList1))
end


