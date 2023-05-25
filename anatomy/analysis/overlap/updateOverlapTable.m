% -------------------------------------------------------------------------
% function to update one row or column of an overlap table with new
% calculations
%
% EXAMPLE USAGE:
%{
overlapPath = ['D:\Fly Imaging\Erica Interneuron Stacks\' ...
    'overlap_calculations\masked\'] ;
overlapFn = 'IN_DN_overlap_table.mat' ;
overlapPathFull = fullfile(overlapPath, overlapFn) ;
overlap_table = importdata(overlapPathFull) ;

neuron_to_update = 'DNp01' ;

maskedFlag = true ;
symFlag = true ;

overlap_table_updated = updateOverlapTable(overlap_table, neuron_to_update, maskedFlag, symFlag) ;

overlap_table = overlap_table_updated ;
save(overlapPathFull, 'overlap_table')
%}
%
% -------------------------------------------------------------------------
function overlap_table = updateOverlapTable(overlap_table, neuron_to_update,...
    maskedFlag, symFlag, dataDim)
% ------------------------------
% inputs/params
if ~exist('maskedFlag','var') || isempty(maskedFlag)
    maskedFlag = true ;
end
if ~exist('symFlag','var') || isempty(symFlag)
    symFlag = true ;
end
if ~exist('dataDim','var') || isempty(dataDim)
    dataDim = [] ;
end

% indices for mask cell array
INPUT = 1 ;
OUTPUT = 2 ;

% data type to load
if symFlag
    dataType = 'coords_sym' ;
    maskType = 'mask_coords_sym' ;
else
    dataType = 'coords_non_sym' ;
    maskType = 'mask_coords_non_sym' ;
end
% -------------------------------------------------------------------
% first determine whether "neuron_to_update" is in the table rows or
% columns
rowNames = overlap_table.Properties.RowNames ;
colNames = overlap_table.Properties.VariableNames ;

row_match_idx = cellfun(@(y) strcmpi(neuron_to_update, y), rowNames) ;
col_match_idx = cellfun(@(y) strcmpi(neuron_to_update, y), colNames) ;

if isempty(dataDim) || (~strcmpi(dataDim,'row') && ~strcmpi(dataDim,'col'))
    if (sum(row_match_idx) ~= 1) && (sum(col_match_idx) ~= 1)
        fprintf('Could not find match for %s \n', neuron_to_update)
        keyboard
    elseif (sum(row_match_idx) == 1) && (sum(col_match_idx) ~= 1)
        dataDim = 'row' ;
    elseif (sum(row_match_idx) ~= 1) && (sum(col_match_idx) == 1)
        dataDim = 'col' ;
    else
        fprintf('Present in both rows and cols... \n')
        keyboard
    end
end
% ----------------------------------------------------------------
% based on whether it's in a row or column, calculate new overlap

% find names of neurons to get overlap against
switch dataDim
    case 'row'
        neuronList = colNames ;
        updateMaskInd = INPUT ;
        otherMaskInd = OUTPUT ;
    case 'col'
        neuronList = rowNames ;
        updateMaskInd = OUTPUT ;
        otherMaskInd = INPUT ; 
    otherwise
        fprintf('Error: did not get data dimension ... \n')
        keyboard
end

% initialize storage array for new overlap calculations
N_cells = length(neuronList) ;
overlap_new = nan(N_cells,1) ;

% if we're dealing with INs, need to account for ID being in table row/col
% label
neuron_to_update_split = strsplit(neuron_to_update,'_') ;
if length(neuron_to_update_split) > 1
    neuron_name = strjoin(neuron_to_update_split(1:end-1),'-') ;
else
    neuron_name = neuron_to_update ;
end

% load data for "neuron_to_update"
neuronToUpdateData = myLoadNeuronData(neuron_name, dataType) ;

% also load/apply mask for neuron to update
if maskedFlag
    % load mask
    neuronToUpdateMask = myLoadNeuronData(neuron_name, maskType) ;
    % get input or output mask (as appropriate)
    neuronToUpdateMask = neuronToUpdateMask{updateMaskInd} ;
    
    % apply mask
    neuronToUpdateData = intersect(neuronToUpdateData,neuronToUpdateMask) ;
end

% ---------------------------------------------------------
% loop over neurons to re-calculate overlap with
neuronType = [] ;
for k = 1:N_cells
    % load data for current neuron
    table_name_curr = neuronList{k} ;
    table_name_split = strsplit(table_name_curr,'_') ;
    if length(table_name_split) > 1
        name_curr = strjoin(table_name_split(1:end-1),'-') ;
    else
        name_curr = table_name_curr ;
    end
    
    [neuronDataCurr, neuronType] = myLoadNeuronData(name_curr, dataType, ...
        neuronType) ;
    
    % apply mask?
    if maskedFlag
        maskDataCurr = myLoadNeuronData(name_curr, maskType, neuronType) ;
        maskDataCurr = maskDataCurr{otherMaskInd} ;

        neuronDataCurr = intersect(neuronDataCurr, maskDataCurr) ;
    end
    
    % get overlap
    overlap_curr = numel(intersect(neuronToUpdateData, neuronDataCurr)) ;
    overlap_new(k) = overlap_curr ;
end

% ------------------------------------------------------------------------
% add new overlap calculations to overlap_table
switch dataDim
    case 'row'
        overlap_table{row_match_idx,:} = overlap_new' ;
        
    case 'col'
        overlap_table{:,col_match_idx} = overlap_new ;
        
    otherwise
        fprintf('Error: did not get data dimension ... \n')
        keyboard
end

end