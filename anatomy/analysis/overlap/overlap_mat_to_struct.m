% -------------------------------------------------------------------------
% script to re-save overlap calculations so that we include labels. this
% should make integrating new images easier
% -------------------------------------------------------------------------
function overlap_struct = overlap_mat_to_struct(overlapPathFull, saveFlag, ...
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
%% fill in structure with data from matrix and labels from directories
overlap_struct = struct() ; 
cc = 1 ; % initialize counter
for i = 1:size(overlap_mat,1)
    for j = 1:size(overlap_mat,2)
    % get current overlap value and check if it's a number
    val_curr = overlap_mat(i, j) ; 
    if isnan(val_curr)
        continue 
    else
        % read labels
        label_1 = dataDir1(i).name ;
        label_2 = dataDir2(j).name ;
        
        % add to structure
        overlap_struct(cc).label_1 = label_1 ;
        overlap_struct(cc).label_2 = label_2 ;
        overlap_struct(cc).overlap_val = val_curr ;
        
        cc = cc + 1 ; % increment counter
    end
    end
end

% -------------------------------------------------------------------------
%% save results?
savePath = fullfile(dataRoot, 'overlap_calculations') ; 
if saveFlag && exist(savePath,'dir')
   saveName = [strjoin(overlap_fn_split(1:(end-1)),'_') '_struct.mat'] ; 
   save(fullfile(savePath, saveName),'overlap_struct')
end

end