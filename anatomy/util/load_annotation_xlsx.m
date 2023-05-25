% -------------------------------------------------------------------------
% function to load annotation data from xlsx file
% -------------------------------------------------------------------------
function [annotationVals, neuron_labels, target_labels, ID_nums, ...
    hemilineage_ID] = load_annotation_xlsx(dataPath, target_labels, ...
    removeDuplicates, neuronNameField, excludeVUMFlag)
% input of which target vnc (or brain, i guess) regions we're looking at
if ~exist('target_labels','var') || isempty(target_labels)
    % list target labels 
    target_labels = {'NTct', 'WTct', 'HTct', 'IntTct', 'LTct', 'T1LegNp', ...
        'T2LegNp','T3LegNp', 'AMNp', 'ANm', 'VAC', 'mVAC'} ; %, 'nonVNC'} ; 
end
% input to determine whether or not we remove neurons with same lineage AND
% innervation pattern (should probably be a yes)
if ~exist('removeDuplicates','var') || isempty(removeDuplicates)
    removeDuplicates = false ;
end
% column name for neuron_labels
if ~exist('neuronNameField','var') || isempty(neuronNameField)
   neuronNameField = 'newName' ;  
end
% whether or not to include VUM
if ~exist('excludeVUMFlag','var') || isempty(excludeVUMFlag)
   excludeVUMFlag = false ;  
end

% -------------------------------------------------------------------------
% load file as table
try
    % detectImportOptions is causing strange errors at the moment...
    opts = detectImportOptions(dataPath);
    annotationsTable = readtable(dataPath,opts) ;
catch
    annotationsTable = readtable(dataPath) ; 
end

% find column indices for targets
target_ind = zeros(1,length(target_labels)) ; 
varNames = annotationsTable.Properties.VariableNames ;
for i = 1:length(target_labels)
    target_ind(i) = find(cellfun(@(y) strcmpi(y, target_labels{i}), ...
        varNames)) ;
    
end
label_ind = find(cellfun(@(y) strcmp(y, neuronNameField), varNames)) ; 

% read in xlsx file info
annotationVals = table2array(annotationsTable(1:end, target_ind)) ;
neuron_labels = annotationsTable.(neuronNameField) ;
ID_nums = annotationsTable.ID ;
hemilineage_ID = annotationsTable.Lineage ; 


% remove duplicate entries?
if removeDuplicates
    % find rows with unique innervation and lineage
   subTable = annotationsTable(1:end,[target_ind, label_ind]) ; 
   [tableUnique, idx_unique] = unique(subTable,'stable') ; 
   
   % take unique annotations and labels
   annotationVals = table2array(tableUnique(1:end,1:length(target_ind))) ; 
   neuron_labels = tableUnique.(neuronNameField) ;
   
   % keep "temporal ID" values only for unique rows
   neuron_labels = neuron_labels(idx_unique) ;
   
   fprintf('Not sure this still works how it should...\n')
   keyboard
end

% exclude empty rows
empty_row_idx = cellfun(@(y) isempty(y), neuron_labels) ; 
neuron_labels = neuron_labels(~empty_row_idx) ; 
annotationVals = annotationVals(~empty_row_idx,:) ; 
ID_nums = ID_nums(~empty_row_idx) ; 
hemilineage_ID = hemilineage_ID(~empty_row_idx) ; 

% -------------------------------
% take out VUM neurons?
if excludeVUMFlag
    % get indices for VUM/PSI neurons
   vum_idx = cellfun(@(y) contains(y, 'VUM','IgnoreCase',true) | ...
       contains(y, 'VPM','IgnoreCase',true), neuron_labels) ;  
   
   % remove these neurons from lists
   annotationVals = annotationVals(~vum_idx,:) ; 
   neuron_labels = neuron_labels(~vum_idx) ; 
   ID_nums = ID_nums(~vum_idx) ; 
   hemilineage_ID = hemilineage_ID(~vum_idx) ; 

end
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% HELPER FUNCTIONS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % ------------------------------------------------------------------------
% %% convert new annotation characters to int
% function vals_out = convert_new_annotation_vals(vals_in, capitalFactor)
% % inputs and params
% if ~exist('capitalFactor','var') || isempty(capitalFactor)
%     % difference between upper and lowercase letters in new annotation
%     % format -- Erica uses e.g. D and d to signify major and minor presence
%     % of dendrites in a region. This factor determines the numerical
%     % difference between "D" and "d"
%    capitalFactor = 2 ;  
% end
% 
% 
% 
% end
