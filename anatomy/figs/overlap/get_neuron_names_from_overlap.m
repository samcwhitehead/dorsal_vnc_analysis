% -------------------------------------------------------------------------
% function to get neuron labels for rows/columns of overlap plot (trying to
% write functions to clean up that plotting code)
%
% "dim" is either 1 or 2 (rows or columns)
% 
% NB: this only works if overlap data is in table format
% -------------------------------------------------------------------------
function labels = get_neuron_names_from_overlap(overlap_data, dim)
% -------------------------------
%% inputs 
% set default dim value
if ~exist('dim','var') || isempty(dim)
    fprintf('No dimension specified -- getting row labels \n')
    dim = 1 ;  
end
% check that overlap data is in correct format
if ~istable(overlap_data)
    fprintf('Error: input data must be a table to get labels from! \n') 
    labels = {} ; 
    return
end

% -----------------------------------------
%% if input data is okay, get labels
% read out table row or column labels
if dim == 1
    tableNames = overlap_data.Properties.RowNames ;
elseif dim == 2
    tableNames = overlap_data.Properties.VariableNames ;
else
    fprintf('Error: invalid table dimension selected: %d \n', dim)
    labels = {} ;
    return
end

% initialize storage for labels
N_names = length(tableNames) ; 
labels = cell(N_names, 1) ;

% loop through and process full names to extract just neuron name 
%  (don't take ID, etc)
for i = 1:N_names
    tableNameCurr = tableNames{i} ;
    tableNameSplit = strsplit(tableNameCurr,'_') ;
    if (length(tableNameSplit) > 1)
        % remaining underscores switched to hyphens
        newName = strjoin(tableNameSplit(1:end-1),'-') ; 
    else
        newName = tableNameCurr ;
    end
    labels{i} = newName ; 
end

end