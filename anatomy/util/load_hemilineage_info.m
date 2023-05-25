% -------------------------------------------------------------------------
% function to load hemilineage table info
% -------------------------------------------------------------------------
function hemilineage_table = load_hemilineage_info(dataPath)
% input(s)
if ~exist('dataPath','var') || isempty(dataPath)
    dataPath = fullfile('D:\Dropbox\Paper Manuscripts', ...
        'Janelia vnc paper', 'cell_lists', 'hemilineage_table.xlsx') ; 
end

% -----------------------------------------------------------
% load data file 
try
    % detectImportOptions is causing strange errors at the moment...
    opts = detectImportOptions(dataPath);
    hemilineage_table = readtable(dataPath,opts) ;
catch
    hemilineage_table = readtable(dataPath) ; 
end 

end
