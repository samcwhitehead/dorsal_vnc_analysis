% -------------------------------------------------------------------------
% function to remove rows/columns from a table (to be used e.g. when a
% neuron is marked as duplicate and removed from the list)
% -------------------------------------------------------------------------
function table_out = removeNeuronFromTable(table_in, neuron_name) 
% initialize table out as copy of table in
table_out = table_in ; 

% find row and/or column (variable) names that match neuron to remove
rowNames = table_out.Properties.RowNames ; 
varNames = table_out.Properties.VariableNames ; 

row_match_idx = cellfun(@(y) strcmp(y, neuron_name), rowNames) ; 
var_match_idx = cellfun(@(y) strcmp(y, neuron_name), varNames) ; 

% indices for rows/columns that we should keep
good_row_idx = ~row_match_idx ; 
good_var_idx = ~var_match_idx ; 

% prune table_out
table_out = table_out(good_row_idx, good_var_idx) ; 

end