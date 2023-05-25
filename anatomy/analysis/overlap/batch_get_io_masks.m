% -------------------------------------------------------------------------
% script to get input/output masks for a set of neurons
% -------------------------------------------------------------------------
%% options
dataTypeList = {'VUM'} ; %{'IN', 'DN', 'MN'} ; % 'IN'

% how should masks be stored?
output_type = 'coords' ; % 'coords' (coordinate indices) | 'bw' (binary images)

% save over?
overWriteFlag = false ;

% print after each neuron is finished?
verboseFlag = true ;

% use new annotations?
newAnnotFlag = true ;
capitalScale = 1 ;  % some params for new annotations
useIpsiContraFlag = false ; 

% what to append to filenames when saving
fileSuffixStr = '' ;
if strcmp(output_type, 'coords')
    fileSuffixStr = [fileSuffixStr '_coords'] ;
end

% where to look for neuropil images (registered/segmented)
[mfilePath, ~, ~] = fileparts(mfilename('fullpath')) ;
figDirectory = fileparts(mfilePath) ;
parentDirectory = fileparts(figDirectory) ;
imPath = fullfile(parentDirectory, 'data', 'processed_anatomy_data') ;
neuropilPath = fullfile(imPath, 'VNC neuropil') ;
neuropil_dir = dir(fullfile(neuropilPath, '*.mat')) ;

% find indices of "coords" filenames
coords_file_idx = arrayfun(@(x) contains(x.name,'coords'), neuropil_dir) ;
if strcmp(output_type,'coords')
    neuropil_dir = neuropil_dir(coords_file_idx) ;
else
    neuropil_dir = neuropil_dir(~coords_file_idx) ;
end

% target neuropil region list
target_labels = {'NTct', 'WTct', 'HTct', 'IntTct1', 'IntTct2','IntTct3',...
    'LTct1', 'LTct2', 'LegNp1', 'LegNp2','LegNp3', 'AMNp', 'ANm'} ; %, 'nonVNC'} ;

% ------------------------------------------------------------------
%% loop over data types
for d = 1:length(dataTypeList)
    
    % ---------------------------------------------------------------------
    % get path info for current data type (NB: all we really need is neuron
    % name for each, so we could use pretty much any directory)
    dataType = dataTypeList{d}  ;
    dataRoot = fullfile(imPath, dataType) ;
    dataPath = fullfile(dataRoot, 'binarized_sym') ;
    
    % where to save output
    savePath = fullfile(dataRoot, 'io_masks') ;
    if ~exist(savePath,'dir')
        mkdir(savePath)
    end
    
    % get directory of bw images for current data type
    bwDir = dir(fullfile(dataPath, '*bw.mat')) ;
    N_bw = length(bwDir) ;
    
    % -------------------------------------------------------
    %% load annotation data for current data type
    % NB: depends on whether or not we're using new or old annotations
    % (except for DNs I guess, which don't have new annotations
    annPath = fullfile(parentDirectory, 'cell_annotations') ;
    switch dataType
        case 'DN'
            annFn = 'DN_annotations.xlsx'  ; % data file containing annotations
        case {'IN', 'VUM', 'MN'}
            annFn = 'IN_annotations.xlsx' ;
        otherwise
            fprintf('Error: unrecognized neuron type: %s \n', dataType)
            keyboard
    end
    
    annotationPathFull = fullfile(annPath, annFn) ;
    [annotationVals, neuron_labels, target_labels, ~] = ...
        load_annotation_xlsx(annotationPathFull, target_labels) ;
    
    % if using new annotations, convert from cell array to 1,2,3 matrix
    if newAnnotFlag && ~strcmp(dataType,'DN')
        % turn cell array into NxMx(2 or 4) array (2 if we're not using
        % ipsi/contra distinction). NB: currently the useIpsiContra option
        % is not supported, since we haven't made sure neurons all have
        % cell bodies on the same side/haven't made ipsi/contra masks
        [~, annotationValsMat] = convertNewAnnotationVals(annotationVals,...
            capitalScale, useIpsiContraFlag) ; 
        
        % get dimensions of converted annotations and initialize output,
        % which will be a 2D array with standard coding of 1 = input, 
        % 2 = output, 3 = both
        mat_size = size(annotationValsMat) ; 
        annotationValsNew = zeros(mat_size(1), mat_size(2)) ; 
        
        % loop over 1st and 2nd dims to convert entries
        for rnum = 1:mat_size(1) 
            for cnum = 1:mat_size(2)
                annot_curr = squeeze(annotationValsMat(rnum, cnum,:)) ; 
                annotationValsNew(rnum, cnum) = ...
                    newAnnotationVecToNum(annot_curr, capitalScale) ;
            end
        end
        
        % redefine annotationVals with new 2D array
        annotationVals = annotationValsNew ; 
    end
    
    % convert neuron labels from field-appropriate form to label form
    neuron_labels = cellfun(@(y) strrep(y,'_','-'), neuron_labels,...
        'UniformOutput', false) ;
    % ----------------------------------------------
    %% loop over images and get masks for them
    for ind = 1:N_bw
        tic
        
        fprintf('Getting masks for %d / %d stacks \n', ind, N_bw)
        % get filename for image
        dataFilename = fullfile(bwDir(ind).folder, bwDir(ind).name) ;
        [~, fn, ext] = fileparts(dataFilename) ;
        % set save path
        maskSavePath = fullfile(savePath, [fn '_mask' fileSuffixStr '.mat']) ;
        
        % check if we've already converted this one
        if (~overWriteFlag) && exist(maskSavePath,'file')
            fprintf('Already analyzed file: %s \n', dataFilename)
            continue
        end
        
        % -----------------------------------------------------------
        % get neuron name from dataFilename
        dataFilenameSplit = strsplit(fn,'_') ;
        neuron_name = dataFilenameSplit{1} ;
        
        % --------------------------------------------
        % calculate masks
        neuropil_masks = getNeuropilMasks(neuron_name, annotationVals,...
            neuron_labels, target_labels, neuropil_dir, output_type, ...
            verboseFlag) ;
        
        % -------------------------------------------
        % save masks
        save(maskSavePath,'neuropil_masks')
        
        toc
    end
end
