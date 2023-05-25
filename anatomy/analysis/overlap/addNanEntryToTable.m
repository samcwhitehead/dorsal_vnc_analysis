% -------------------------------------------------------------------------
% function to add an empty row/col to overlap table (overlap_data)
% -------------------------------------------------------------------------
function overlap_data = addNanEntryToTable(overlap_data, neuron_name, dim)
% get table size
[nrows, ncols] = size(overlap_data) ; 

% also read out row and column names (from unaltered table)
rowNames = overlap_data.Properties.RowNames ; 
colNames = overlap_data.Properties.VariableNames ; 

% ---------------------------------------------------
% add an empty placeholder row/col for new entry
switch dim
    case 1
        % if we're adding a row, create a row vector table of nan values 
        % to add to overlap_data
        new_entry = nan(1,ncols) ; 
        new_entry = array2table(new_entry, 'VariableNames', colNames,...
            'RowNames',{neuron_name}) ; 
        overlap_data = [overlap_data ; new_entry] ; 
        
        % sort overlap_data rows so that new entry is in the right order
        rowNames = overlap_data.Properties.RowNames ; % NEW row names
        [~, sort_ind] = sort(rowNames) ; 
        overlap_data = overlap_data(sort_ind,:) ; 
        
    case 2
        % if we're adding a column, create a column vector table of nan
        % values to add to overlap_data
        new_entry = nan(nrows,1) ; 
        new_entry = array2table(new_entry,'VariableNames', {neuron_name},...
            'RowNames', rowNames) ; 
        overlap_data = [overlap_data, new_entry] ; 
        
        % sort overlap_data cols so that new entry is in the right order
        colNames = overlap_data.Properties.VariableNames ;  % NEW col names
        [~, sort_ind] = sort(colNames) ; 
        overlap_data = overlap_data(:, sort_ind) ; 
    otherwise
        fprintf('Invalid dim entry: %d \n', dim)
        keyboard
end

end