% -------------------------------------------------------------------------
% function to pull out shortened name for each stack file that can be used
% as unique identifier in tables, etc
% -------------------------------------------------------------------------
function [lineName, justName, justID] = getLineName(fileName, dataType,...
    forTableFlag, wrongFlag)
% -----------------------------------------
%% inputs
if ~exist('forTableFlag','var') || isempty(forTableFlag)
   forTableFlag = false ; % use this to replace any dashes with underscores (dashes are invalid in table var names) 
end
if ~exist('wrongFlag','var') || isempty(wrongFlag)
   wrongFlag = false ; % get the wrong name? (to be used for renaming)
end
% -----------------------------------------
%% search strings (depending on data type)
% assume we won't search ID number, but change if we have IN or VUM
searchIdFlag = false ;
% use data type to determine what to search for in filename
switch dataType
    case 'IN'
        if wrongFlag
            line_search_str = '(?<lineName>[NWHXILBUMO]{1,4}([0-9]|-){0,5}_)' ; 
        else
            line_search_str = '(?<lineName>[NWHXILBUMOT2PS]{1,4}([0-9]|-){0,7}_)' ; %'(?<lineName>[NWHVUMOIPSMXBL]{1,4}([0-9]|-){0,5}_)' ;
        end
        id_search_str = '(?<lineID>_[^_]*?\d{4}[a-zA-Z]{0,5})' ; 
        searchIdFlag = true ;
    case 'VUM'
        line_search_str = '(?<lineName>([TVUMAPabcde1-9]{1,6}|PSI)_)' ;
        id_search_str = '(?<lineID>_[^_]*?(\d{4}[a-zA-Z]{0,5})|V\d{2,3})' ; 
        searchIdFlag = true ;
    case 'DN' 
        line_search_str = '(?<lineName>DN\w{1,6}_)' ; 
    case 'MN'
        line_search_str = '(?<lineName>[DVLMshgitpNnb]{1,4}\w+_)' ; 
    otherwise
        fprintf('Invalid dataType: %s \n', dataType) 
        keyboard
end

% ----------------------------------------------------
%% search for expressions in string
lineSearchResults = regexp(fileName, line_search_str,'names') ;
lineName = lineSearchResults(1).lineName ;

% remove extra underscores
lineName = erase(lineName,'_') ;

if forTableFlag
   lineName = strrep(lineName,'-','_') ;  
end

% store plain name (without ID)
justName = lineName ; 

% search for line ID?
if searchIdFlag
    iDsearchResults =  regexp(fileName, id_search_str,'names') ;
    lineID = iDsearchResults(1).lineID ;
    % remove extra underscores
    lineID = erase(lineID,'_') ;
    
    % add ID to line name
    lineName = strjoin({lineName, lineID},'_') ;
    
    % store plain line ID
    justID = lineID ; 
else
    justID = '' ; 
end

end