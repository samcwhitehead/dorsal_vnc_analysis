% -------------------------------------------------------------------------
% script to re-save overlap calculations so that we include labels. this
% should make integrating new images easier
%
%{
overlapPathFull = fullfile(['D:\Fly Imaging\Erica Interneuron ' ...
    'Stacks\overlap_calculations\'], 'MN_MN_overlap_mat.mat')  ;
overlap_table = overlap_mat_to_table(overlapPathFull) ; 
%}
% -------------------------------------------------------------------------
function overlap_table = overlap_mat_to_table(overlapPathFull, saveFlag, ...
    dataRoot, bwFolderName) 
% ---------------------
%% inputs and params
if ~exist('saveFlag','var') || isempty(saveFlag)
   saveFlag = true ;  % save output structure
end
if ~exist('dataRoot','var') || isempty(dataRoot)
   dataRoot = 'D:\Fly Imaging\Erica Interneuron Stacks\' ;  % path to image data
end
if ~exist('bwFolderName','var') || isempty(bwFolderName)
   bwFolderName = 'binarized_new' ;  % name of folder where bw mat files are stored
end
% -----------------------------------------------------------------------
%% get info on sets of neurons being compared based on overlap file path
[~, overlap_fn, ~] = fileparts(overlapPathFull) ; 
overlap_fn_split = strsplit(overlap_fn,'_') ; 

dataName1 = overlap_fn_split{1} ; 
dataName2 = overlap_fn_split{2} ; 
if ~strcmp(overlap_fn_split{3}, 'overlap')
   fprintf('Error determining group names \n')
   keyboard
end

dataPath1 = fullfile(dataRoot, dataName1,bwFolderName) ; 
dataPath2 = fullfile(dataRoot, dataName2,bwFolderName) ; 

% -------------------------------------------------------------------------
%% load overlap mat and get names of files in the two data directories
overlap_mat = importdata(overlapPathFull) ; 
dataDir1 = dir(fullfile(dataPath1,'*bw.mat')) ; 
dataDir2 = dir(fullfile(dataPath2,'*bw.mat')) ; 

% get indices of non-nan entries in overlap matrix
%[rows, cols] = find(~isnan(overlap_mat)) ; 

% -------------------------------------------------------------------------
%% fill in table with data from matrix and labels from directories
overlap_table = array2table(overlap_mat) ; 

labelList1 = cell(size(overlap_mat,1),1) ; 
for i = 1:size(overlap_mat,1)
    labelList1{i} = getLineName(dataDir1(i).name, dataName1, true) ; 
end
labelList2 = cell(size(overlap_mat,2),1) ; 
for j = 1:size(overlap_mat,2)
    labelList2{j} = getLineName(dataDir2(j).name, dataName2, true) ; 
end
overlap_table.Properties.RowNames = labelList1 ; 
overlap_table.Properties.VariableNames = labelList2 ; 


% -------------------------------------------------------------------------
%% save results?
savePath = fullfile(dataRoot, 'overlap_calculations') ; 
if saveFlag && exist(savePath,'dir')
   saveName = [strjoin(overlap_fn_split(1:(end-1)),'_') '_table.mat'] ; 
   save(fullfile(savePath, saveName),'overlap_table')
end

end