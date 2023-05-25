% -------------------------------------------------------------------------
% script to run clustering on manual annotations
% -------------------------------------------------------------------------
% using DNs or INs+VUMs (default)?
DNFlag = false ; 

% path to current script, used to locate data and save folder
[mfilePath, ~, ~] = fileparts(mfilename('fullpath')) ; 
analysisPath = fileparts(mfilePath) ; 
parentDirectory = fileparts(analysisPath) ;
rootPath = fullfile(parentDirectory, 'data', 'cell_annotations') ;
savePath = fullfile(parentDirectory,'data','cluster_output') ;

% flag options
saveFlag = false ; 
overWriteFlag = true ; 
debugFlag = true ; 
excludeVUMFlag = false ; % new annotation file does not contain VUMs 

% make sire save directory exists
if saveFlag && ~exist(savePath,'dir')
    mkdir(savePath)   
end

% different case depending on whether we're using DNs or INS+VUMs
if DNFlag
    % data file containing annotations
    dataFilename = 'DN_annotations.xlsx'  ; 
else
    dataFilename = 'IN_annotations.xlsx'  ; 
end

% target labels in annotation files
target_labels = {'NTct', 'WTct', 'HTct', 'IntTct1', 'IntTct2','IntTct3',...
    'LTct1', 'LTct2', 'LegNp1', 'LegNp2','LegNp3', 'AMNp', 'ANm', 'Y'} ; %, 'nonVNC'} ;

% clustering options
if DNFlag
    max_k = 3 ;
else
    max_k = 29 ;  %21 ; % 16 
end
distanceType = 'correlation' ; 
linkageMethod = 'complete' ; 

% parameters for new annotations
capitalScale = 2 ; % how much to weight upper vs lowercase letters
useIpsiContraFlag = true ; % preserve information on ipsi vs contra annotations
rmWeakAnnotFlag = true ; % only pay attention to capital letters 

% ---------------------------------------------------------------------
% load annotation data
dataPath = fullfile(rootPath, dataFilename) ;
[annotationVals, neuron_labels, target_labels, ID_nums, hemilineage_IDs] = ...
    load_annotation_xlsx(dataPath,target_labels,[],[],excludeVUMFlag) ;

% ---------------------------------------------------------------------
% run clustering function
[conn_clust_struct, conn_clust_table] = ...
    cluster_io_annotations(annotationVals, neuron_labels, ID_nums, ...
    target_labels, max_k, distanceType, linkageMethod, debugFlag, ...
    capitalScale, useIpsiContraFlag, rmWeakAnnotFlag) ;

% ---------------------------------------------------------------------
% add hemilineage ID to clustering struct
for k = 1:length(conn_clust_struct)
    % get current neuron ID (if it's not empty)
    ID_curr = conn_clust_struct(k).ID ; 
    if isempty(ID_curr)
       continue 
    end
    % find matching index for data loaded from xlsx file
    match_idx = cellfun(@(y) strcmp(y, ID_curr), ID_nums) ; 
    
    if sum(match_idx) ~= 1
        fprintf('Error: could not find match for neuron ID %s \n', ID_curr)
        keyboard
    end
    
    % assign hemilineage ID to the matching element in conn_clust_struct
    conn_clust_struct(k).hemilineageID = hemilineage_IDs{match_idx} ;
end

% ---------------------------------------------------------------------
% store original cell array output for convenience
conn_clust_struct(1).annotationValsCell = annotationVals ;  

% ---------------------------------------------------------------------
% save output?
if saveFlag
    % save path 
    saveFn = 'conn_clust_struct' ; 
    suffixStr = '' ; 
    if DNFlag 
       suffixStr = [suffixStr '_DN'] ;
    end
   
    suffixStr = sprintf('%s_capScale_%d',suffixStr,capitalScale) ;
    if useIpsiContraFlag
        suffixStr = sprintf('%s_useIpsiContra',suffixStr) ;
    end
    if rmWeakAnnotFlag
        suffixStr = sprintf('%s_rmWeak',suffixStr) ;
    end
  
    savePathFull = fullfile(savePath,sprintf('%s%s.mat',saveFn,suffixStr)) ; 
    
    % save struct?
    if ~exist(savePathFull,'file') || overWriteFlag
        save(savePathFull, 'conn_clust_struct')
        % save table
        writetable(conn_clust_table, ...
             fullfile(savePath,sprintf('%s%s.xlsx',saveFn,suffixStr)))
    end
    
    if debugFlag
       debug_fig_fn = sprintf('cluster_criteria_opt%s',suffixStr) ; 
       savefig(gcf, fullfile(savePath, [debug_fig_fn '.fig']))
       exportgraphics(gcf, fullfile(savePath, [debug_fig_fn '.png']),...
           'Resolution',300)
    end
end