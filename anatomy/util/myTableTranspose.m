% -------------------------------------------------------------------------
% quick function to take table transpose -- i'm sure there's an easier way
% to do this, but whatever
% -------------------------------------------------------------------------
function table_data_transpose = myTableTranspose(table_data)
% i guess first check if we indeed have a table?
if ~istable(table_data)
    fprintf('Error: input is not a table \n')
    try 
        table_data_transpose = table_data' ; 
    catch
        fprintf('Could not take transpose \n') 
        table_data_transpose = [] ; 
    end
    return
end

% then use the built-in function...
table_data_transpose = rows2vars(table_data) ;

% ...but need to deal with the row and column labels
originalVariableNames = table_data_transpose.OriginalVariableNames ;
table_data_transpose = table_data_transpose(:,2:end) ;
table_data_transpose.Properties.RowNames = originalVariableNames ;
        
end