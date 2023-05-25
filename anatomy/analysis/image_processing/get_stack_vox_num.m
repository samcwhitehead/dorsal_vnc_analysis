% -------------------------------------------------------------------------
% quick and dirty script to get voxel count for each stack
% -------------------------------------------------------------------------
dataRoot = 'D:\Fly Imaging\Erica Interneuron Stacks\VUM\' ;
coordsDir = dir(fullfile(dataRoot,'binarized_new', '*_coords.mat')) ;

saveFlag = true ;
maskedFlag = true ; 
N_stacks = length(coordsDir) ; 

% determine which data type we're looking at by folder name
dataRootSplit = strsplit(dataRoot,'\') ; 
empty_idx = cellfun(@(y) isempty(y), dataRootSplit) ; 
dataRootSplit = dataRootSplit(~empty_idx) ; 
dataType = dataRootSplit{end} ; 

% ------------------------------
% initialize output
if maskedFlag
    vox_size_vec = nan(N_stacks,2) ; 
    maskDir = dir(fullfile(dataRoot, 'io_masks', '*_mask_coords.mat')) ; 
    INPUT = 1 ; 
    OUTPUT = 2 ; 
else
    vox_size_vec = nan(N_stacks,1) ; 
end
lineNameList = cell(N_stacks,1) ; % list of line names corresponding to vox nums

% --------------------------------------
% loop through bw images and get size
for i = 1:N_stacks
    tic
    
    % ----------------------------------------------------------
    % first get name identifiers for this line
    fn = coordsDir(i).name ; 
    lineName = getLineName(fn, dataType) ; 
    
    % ------------------------------------------------------
    % load bw data and count voxels
    fn_full = fullfile(coordsDir(i).folder, coordsDir(i).name) ; 
    coords_curr = importdata(fn_full) ; 
    
    % ---------------------------------------------
    % store results
    if maskedFlag
        % if using neuropil region masks, load current mask
        mask_fn = fullfile(maskDir(i).folder, maskDir(i).name) ; 
        mask_curr = importdata(mask_fn) ; 
        
        % store input and output total voxel count
        vox_size_vec(i,INPUT) = numel(intersect(mask_curr{INPUT}, ...
            coords_curr)) ;
        vox_size_vec(i,OUTPUT) = numel(intersect(mask_curr{OUTPUT}, ...
            coords_curr)) ;
    else 
        vox_size_vec(i) = numel(coords_curr) ; 
    end
    
    % either way, store line name
    lineNameList{i} = lineName ; 
    
    toc
end

% -------------------------------
% generate table from output
vox_size_table = table(lineNameList,vox_size_vec) ; 
% rename columns
vox_size_table.Properties.VariableNames{1} = 'lineName' ; 
vox_size_table.Properties.VariableNames{2} = 'voxCount' ; 

% ------------------------------
% save results?
if saveFlag
    if maskedFlag
        suffixStr = '_masked' ; 
    else
        suffixStr = '' ; 
    end
   save(fullfile(dataRoot,'binarized_new', ...
       ['vox_size_vec' suffixStr '.mat']), ...
       'vox_size_vec')
   save(fullfile(dataRoot,'binarized_new',...
       ['vox_size_table' suffixStr '.mat']), ...
       'vox_size_table')
end