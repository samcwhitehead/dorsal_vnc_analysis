% -------------------------------------------------------------------------
% script to add or update overlap calculations for a set of neurons to all 
% overlap tables containing this neuron type
% -------------------------------------------------------------------------
% ------------------------------------
%% path and params
% path info
[mfilePath, ~, ~] = fileparts(mfilename('fullpath')) ;
figDirectory = fileparts(mfilePath) ;
parentDirectory = fileparts(figDirectory) ;
rootPath = fullfile(parentDirectory, 'data') ;
    
overlapPath = fullfile(rootPath, 'overlap_calculations') ; 

% which neurons to add, and what type of neuron each is
new_neuron_list = {} ;  % {'T2VUMb_3083'} ; 
new_neurons_data_types = {} ;  % {'VUM'} ; 

% boolean options for "addNeuronsToTable.m"
saveFlag = true ;  % save output?
maskedFlag = true ;  % calculate masked overlap?
symFlag = true ;    % use symmetrized neuron stacks?
% if neurons we're adding already exist in table, skip them?
% (don't need to try this, because we can just update)
overWriteFlag = true ;

% data dim names
data_dim_names = {'row','col'} ; 

% -------------------------------------------------------
%% get overlap directory
% update path to overlap data if using masked data
if maskedFlag
    overlapPath = fullfile(overlapPath, 'masked') ;
end

% get directory of overlap matrices
overlapDir = dir(fullfile(overlapPath, '*_overlap_table.mat')) ; 

% find first and second data types for each overlap data table
overlapFnSplit = arrayfun(@(x) strsplit(x.name,'_'), overlapDir, ...
    'UniformOutput',0) ; 
dataTypeList1 = cellfun(@(y) y{1}, overlapFnSplit, 'UniformOutput',0) ; 
dataTypeList2 = cellfun(@(y) y{2}, overlapFnSplit, 'UniformOutput',0) ; 

% ------------------------------------------------------------
%% loop over new neurons to add
for k = 1:length(new_neuron_list)
   % current neuron and its data type
   neuron_curr = new_neuron_list(k) ;  % keep as cell!
   data_type_curr = new_neurons_data_types{k} ; 
   
   % -----------------------------------------------------------------
   %% find all overlap tables with current neuron type as ROWS
   % ...then add neuron to these tables
   table_row_ind = find(cellfun(@(y) strcmp(data_type_curr, y), ...
       dataTypeList1)) ; 
   neuron_dim = 1 ; 
   % loop over overlap tables with current neuron type as rows
   for m = 1:length(table_row_ind)
       % index and data types
       ind = table_row_ind(m) ; 
       dataType1 = dataTypeList1{ind} ; 
       dataType2 = dataTypeList2{ind} ; 
       
       % probably dumb, but double check that data types match
       if ~strcmp(dataType1, data_type_curr)
          keyboard 
       end
       
       % --------------------------------------------------------------
       % load overlap data table to see if neuron is already present
       overlap_fn = sprintf('%s_%s_overlap_table.mat', dataType1,dataType2);
       overlap_data = importdata(fullfile(overlapPath, overlap_fn)) ; 
       
       % check rows for current neuron
       rowNames = overlap_data.Properties.RowNames ; 
       match_idx = cellfun(@(y) strcmp(neuron_curr{1}, y), rowNames) ; 
       
       % -----------------------------------------------------------------
       % now update or add, depending on whether we found a match
       if sum(match_idx) == 1
           % if we find a match, update overlap table
           fprintf('Updating %s overlap for %s vs %s table ... \n', ...
               neuron_curr{1},  dataType1, dataType2)
           
            overlap_data = updateOverlapTable(overlap_data, neuron_curr{1},...
                maskedFlag, symFlag, data_dim_names{neuron_dim}) ;
           
            % save updated overlap table
            if saveFlag
                save(fullfile(overlapPath, overlap_fn), 'overlap_data') 
            end
            
           fprintf('Completed updating %s overlap for %s vs %s table\n', ...
               neuron_curr{1}, dataType1, dataType2)
       
       elseif sum(match_idx) < 1
           % if no match, add neuron
           fprintf('Adding %s to %s vs %s overlap ... \n', ...
               neuron_curr{1},  dataType1, dataType2)
           
           overlap_data = addNeuronsToTable(neuron_curr, neuron_dim, ...
               dataType1, dataType2, overlapPath, saveFlag, ...
               maskedFlag, symFlag, overWriteFlag) ;
           
           fprintf('Completed adding %s to %s vs %s overlap\n', ...
               neuron_curr{1}, dataType1, dataType2)
       else
           % if multiple matches, then this is probably and error
           fprintf('Error: multiple matches for neuron %s \n', ...
               neuron_curr{1}) ;
           keyboard
       end
   end
   
   % -----------------------------------------------------------------
   %% find all overlap tables with current neuron type as COLUMNS
   % ...then add neuron to these tables
   table_col_ind = find(cellfun(@(y) strcmp(data_type_curr, y), ...
       dataTypeList2)) ; 
   neuron_dim = 2 ; 
   % loop over overlap tables with current neuron type as rows
   for n = 1:length(table_col_ind)
       % index and data types
       ind = table_col_ind(n) ; 
       dataType1 = dataTypeList1{ind} ; 
       dataType2 = dataTypeList2{ind} ; 
       
       % probably dumb, but double check that data types match
       if ~strcmp(dataType2, data_type_curr)
          keyboard 
       end
       
        % --------------------------------------------------------------
       % load overlap data table to see if neuron is already present
       overlap_fn = sprintf('%s_%s_overlap_table.mat', dataType1,dataType2);
       overlap_data = importdata(fullfile(overlapPath, overlap_fn)) ; 
       
       % check COLUMNS for current neuron
       colNames = overlap_data.Properties.VariableNames ; 
       match_idx = cellfun(@(y) strcmp(neuron_curr{1}, y), colNames) ; 
       
       % -----------------------------------------------------------------
       % now update or add, depending on whether we found a match
       if sum(match_idx) == 1
           % if we find a match, update overlap table
           fprintf('Updating %s overlap for %s vs %s table ... \n', ...
               neuron_curr{1},  dataType1, dataType2)
           
           overlap_data = updateOverlapTable(overlap_data, neuron_curr{1},...
                maskedFlag, symFlag, data_dim_names{neuron_dim}) ;
           
             % save updated overlap table
            if saveFlag
                save(fullfile(overlapPath, overlap_fn), 'overlap_data') 
            end
            
           fprintf('Completed updating %s overlap for %s vs %s table\n', ...
               neuron_curr{1}, dataType1, dataType2)
       
       elseif sum(match_idx) < 1
           % if no match, add neuron
           fprintf('Adding %s to %s vs %s overlap ... \n', ...
               neuron_curr{1},  dataType1, dataType2)
           
           overlap_data = addNeuronsToTable(neuron_curr, neuron_dim, ...
               dataType1, dataType2, overlapPath, saveFlag, ...
               maskedFlag, symFlag, overWriteFlag) ;
           
           fprintf('Completed adding %s to %s vs %s overlap\n', ...
               neuron_curr{1}, dataType1, dataType2)
       else
           % if multiple matches, then this is probably and error
           fprintf('Error: multiple matches for neuron %s \n', ...
               neuron_curr{1}) ;
           keyboard
       end
   end
    
end