% -------------------------------------------------------------------------
% FUNCTION to add new neuron(s) to overlap table  
% 
% Example usage:
%{
new_neuron_list = {'hg3'} ; 
new_neuron_dims = 1.*ones(size(new_neuron_list)) ; 

dataType1 = 'MN' ; 
dataType2 = 'MN' ; 

overlap_data = addNeuronsToTable(new_neuron_list, dataType1, ...
    dataType2, dataRoot, saveFlag, maskedFlag, symFlag, overWriteFlag) ; 
%}
%
% -------------------------------------------------------------------------
function overlap_data = addNeuronsToTable(new_neuron_list, ...
    new_neuron_dims, dataType1, dataType2, overlapPath, saveFlag, ...
    maskedFlag, symFlag, overWriteFlag)
% ------------------------------------
%% inputs and params
if ~exist('overlapPath','var') || isempty(overlapPath)
    % path to overlap data storage
    overlapPath = ['D:\Fly Imaging\Erica Interneuron Stacks\'...
    'overlap_calculations\'] ; 
end
if ~exist('saveFlag','var') || isempty(saveFlag)
    saveFlag = false ;  % save output?
end
if ~exist('saveFlag','var') || isempty(saveFlag)
    maskedFlag = true ;  % calculate masked overlap?
end
if ~exist('saveFlag','var') || isempty(saveFlag)
    symFlag = true ;    % use symmetrized neuron stacks?
end
if ~exist('saveFlag','var') || isempty(saveFlag)
    overWriteFlag = false ;  % if neurons we're adding already exist in table, skip them?
end

% corresponding names for neuron overlap matrix dimension values
dim_names = {'row', 'col'} ; % this doesn't need to be changed -- just giving meaning to 1 and 2 above

% --------------------------------------------------------
%% load current overlap data
% path to overlap data -- assumes a fixed location relative to dataRoot
overlapFn = [dataType1 '_' dataType2 '_overlap_table.mat'] ; 
if maskedFlag && ~contains(overlapPath,'masked')
    overlapPath = fullfile(overlapPath, 'masked') ;
end
overlapPathFull = fullfile(overlapPath, overlapFn) ;

% check that overlap data file exists. if so, load it
if ~exist(overlapPathFull, 'file')
    fprintf('Error: could not find file: \n %s \n', overlapPathFull) 
    keyboard
end
overlap_data = importdata(overlapPathFull) ; 

% ---------------------------------------------------------------------
%% loop over neurons to add
for k = 1:length(new_neuron_list)
    % read out name and dimension for current neuron
    neuron_name = new_neuron_list{k} ;
    neuron_dim = new_neuron_dims(k) ; 
    
    % check to see if we actually need to add this neuron
    if neuron_dim == 1
        namesToCheck = overlap_data.Properties.RowNames ;
    elseif neuron_dim == 2
        namesToCheck = overlap_data.Properties.VariableNames ;
    end
    if ismember(neuron_name, namesToCheck) && ~overWriteFlag
       fprintf('Already added %s -- skipping \n', neuron_name) 
       continue 
    end
    
    % add nan row or col to table for new neuron
    overlap_data = addNanEntryToTable(overlap_data,neuron_name,neuron_dim); 
    
    % fill in blank entries with newly calculated overlap
    overlap_data = updateOverlapTable(overlap_data, neuron_name, ...
        maskedFlag, symFlag, dim_names{neuron_dim}) ;

    % print progress
    fprintf('Successfully added %s (%d/%d) to overlap data \n', ...
        neuron_name, k, length(new_neuron_list))
end

% ----------------------------------
%% save results?
if saveFlag
    save(overlapPathFull, 'overlap_data') 
end
end