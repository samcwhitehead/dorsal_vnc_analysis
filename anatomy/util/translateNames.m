% -------------------------------------------------------------------------
%  UNDER CONSTRUCTION
% Function to pull out relevant information from a string describing a
% hemilineage and use it to create plot labels.
%
% Outputs all strings/characters (can convert to numbers later if
% necessary)
% -------------------------------------------------------------------------
function [new_name_str, new_name_cell, hl_num, hl_half, greek_letter, vnc_seg, ...
    neuron_num, driver_name, sex] = translateNames(hl_str, my_delim) 
% ------------------
% params
if ~exist('my_delim','var') || isempty(my_delim)
    my_delim = ' ' ; 
end

%greek_letters = {char(945), char(946), char(947), char(948), char(949) } ; 
greek_letter_list = {'alpha','beta','gamma','delta','epsilon'} ; 
mn_list = {'DLM', 'DVM', 'hg1', 'hg2', 'i1', 'i2', 'ps1', 'tp1', 'tp2', 'tpN'} ; 

translation_struct = struct() ; 
translation_struct.a = char(945) ;
translation_struct.b = char(946) ;
%translation_struct.c = char(947) ;
translation_struct.g = char(947) ;
translation_struct.d = char(948) ;
translation_struct.e =  char(949) ;

% -----------------------------------------------------
% try to determine if there's a delimiter involved
if contains(hl_str, '_') 
    delim = '_' ; 
elseif contains(hl_str, ' ') 
    delim = ' ' ;
else
    delim = '' ; 
end

% -----------------------------------------------------
% keep a running list of string parts we've identified
translated_idx = [] ; 

% -----------------------------------------------------------------
% look for hemilineage number and half. should be in the form XXY
% where XX is an integer and Y is an upper case letter (A or B)
[start_idx, end_idx] = regexp(hl_str, '\d\d[AB]') ;
if isempty(start_idx)
    [start_idx, end_idx] = regexp(hl_str, '\d[AB]') ;
end

% I think hemilineage 17 doesn't have an AB half, so look for that special
% case

if isempty(start_idx) && ~contains(hl_str,'17')
    hl_num = [] ; 
    hl_half = [] ;
    nameCase = 'notHemilineage' ; 
elseif isempty(start_idx) && contains(hl_str,'17')
    hl_num = '17' ; 
    hl_half = [] ;
    nameCase = 'hemilineage' ; 
elseif (numel(start_idx)==1) && isnan(str2double(hl_str(end_idx)))
    hl_num = hl_str(start_idx:(end_idx-1)) ; 
    hl_half = hl_str(end_idx) ; 
    nameCase = 'hemilineage' ; 
    translated_idx = [translated_idx, start_idx:end_idx] ; 
elseif ~isempty(start_idx) && ~isnan(str2double(hl_str(end_idx))) && ...
        ((end_idx - start_idx) == 1)
    hl_num = hl_str(start_idx:end_idx) ; 
    hl_half = '' ; 
    nameCase = 'hemilineage' ;
    translated_idx = [translated_idx, start_idx:end_idx] ; 
else
    disp('problem with hemilineage half')
    keyboard
end

seventeen_idx = strfind(hl_str, '17') ; 
if ~isempty(seventeen_idx)
    translated_idx = [translated_idx, (seventeen_idx):(seventeen_idx + 1)] ; 
end
% --------------------------------------------------------------------
% depending on whether or not this string actually confroms to the
% hemilineage format we expect, do different stuff
switch nameCase
    case 'hemilineage'
        % get greek letter
        [greek_letter, greek_letter_idx] = get_greek_letter(hl_str, ...
            greek_letter_list, translation_struct, translated_idx) ; 
        translated_idx = sort(unique([translated_idx, greek_letter_idx])) ; 

        % get vnc segment
        [vnc_seg, vnc_idx] = get_vnc_segment(hl_str) ; 
        translated_idx = sort(unique([translated_idx, vnc_idx])) ; 

        % get neuron number
        neuron_num = get_neuron_num(hl_str(~translated_idx), delim) ; 

        % get driver name
        driver_name = get_driver_name(hl_str) ; 

        % get fly sex
        sex = get_fly_sex(hl_str) ;
        
        % combine info
        if ~isempty(hl_half) && isempty(hl_num)
            new_name_cell = {hl_num, greek_letter, vnc_seg, ...
                neuron_num, driver_name, sex} ;
        else
            new_name_cell = {[hl_num, hl_half], greek_letter, vnc_seg, ...
                neuron_num, driver_name, sex} ;
        end
        
        % find entries that we couldn't get from string
        empty_idx = cellfun(@(y) isempty(y), new_name_cell) ;
        in_name_idx = false(1,length(new_name_cell)) ;
        in_name_idx([1:4, 6]) = true ;
        
        if (sum(in_name_idx & ~empty_idx) < 1)
            % if this happens, we failed :(
            new_name_str = hl_str ;
        else
            new_name_str = strjoin(new_name_cell(in_name_idx & ~empty_idx), ...
                my_delim) ;
        end
        
    case 'notHemilineage' 
        
        % give empty outputs
        hl_num = [] ;
        hl_half = [] ;
        greek_letter = [] ;
        vnc_seg = [] ;
        neuron_num = [] ;
        driver_name = [] ;
        sex = [] ;
        
        % check if it's a mesVUM line
        if contains(hl_str,'mesVUM','IgnoreCase',true)
            % mesVUM lines are a special case--deal with them
            new_name_cell = cell(1,6) ;
            hl_str_split = strsplit(hl_str,delim) ;
            new_name_str = hl_str_split{1} ;
            return
        end
        % check if it's an MN line
        mn_check = cellfun(@(y) contains(hl_str, y), mn_list) ;
        if (sum(mn_check) > 0)
            new_name_cell = cell(1,6) ;
            new_name_str = hl_str ;  
            return
        end
        % check if it's a "deleted" line
        if contains(hl_str,'deleted','IgnoreCase',true)
            % mesVUM lines are a special case--deal with them
            new_name_cell = cell(1,6) ;
            new_name_str = hl_str ;
            return
        end
        
        % in case none of the above are true, give empty outputs
        new_name_cell = cell(1,6) ;
        new_name_str = [] ;
        
        
    otherwise
        keyboard
end



end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% helper functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------------------------------------------------------
%% search for presence of greek letter. 
% might be abbreviated to just lower case letter, so we'll search for that too
function [greek_letter, greek_letter_idx] = get_greek_letter(hl_str, ...
    greek_letter_list, translation_struct, translated_idx)

% look for greek letters in different forms. could be 'a' 'alpha' or '?'
greek_abbrev_idx = regexp(hl_str, '[a-g]') ;
greek_full_check = false(length(greek_letter_list),1) ; 
greek_char_check = false(length(greek_letter_list),1) ; 
for gl = 1:length(greek_letter_list)
    greek_full_check(gl) = contains(hl_str, greek_letter_list{gl}) ; 
end
translation_fields = fields(translation_struct) ; 
for gc = 1:length(translation_fields)
    greek_char_check(gc) = contains(hl_str, ...
        translation_struct.(translation_fields{gc})) ; 
end

% if we found just character or full name, get indices for where they occur
% in string
if (sum(greek_full_check) == 1)
    greek_letter_idx = strfind(hl_str, greek_letter_list(greek_full_check)) ;
    greek_letter = translation_struct.(greek_letter_list{greek_full_check}(1)) ; 
elseif (sum(greek_char_check) == 1)
    greek_letter_idx = strfind(hl_str, ...
        translation_struct.(translation_fields{greek_char_check})) ; 
    greek_letter = translation_struct.(translation_fields{greek_char_check}) ; 
elseif (numel(greek_abbrev_idx) == 1)
    greek_letter = translation_struct.(hl_str(greek_abbrev_idx)) ; 
    greek_letter_idx = greek_abbrev_idx ; 
elseif (sum(greek_full_check) < 1) && (sum(greek_char_check) < 1) && ...
        (numel(greek_abbrev_idx) < 1)
    greek_letter = [] ; 
    greek_letter_idx = [] ; 
elseif (numel(greek_abbrev_idx) > 1)
    [str_dist, min_idx] = min(abs(translated_idx(end) - greek_abbrev_idx)) ; 
    if (str_dist < 2)
        greek_abbrev_idx = greek_abbrev_idx(min_idx) ;
        greek_letter = translation_struct.(hl_str(greek_abbrev_idx)) ;
        greek_letter_idx = greek_abbrev_idx ;
    else
        greek_letter = [] ;
        greek_letter_idx = [] ;
    end
else
    disp('unsure of what to do about greek letter')
    keyboard
end
% if isempty(greek_abbrev_idx) && (sum(greek_full_check) < 1)
%     greek_letter = '' ; 
%     greek_letter_idx = [] ;
% elseif (sum(greek_full_check) == 1)
%     greek_letter = greek_letter_list{greek_full_check}(1) ; 
% elseif (numel(greek_abbrev_idx) == 1) && (sum(greek_full_check) ~= 1)
%     greek_letter = hl_str(greek_abbrev_idx) ; 
% else
%     disp('unsure of what to do about greek letter')
%     keyboard
% end
% 
% if length(greek_letter) == 1 
%    greek_letter = translation_struct.(greek_letter) ;  
% end
end

% -------------------------------------------------------------------------
%% look for vnc segment
function [vnc_seg, vnc_idx] = get_vnc_segment(hl_str) 
[vnc_start_idx, vnc_end_idx] = regexp(hl_str, '[t]\d') ;
if ~isempty(vnc_start_idx)
    vnc_seg = hl_str(vnc_start_idx:vnc_end_idx) ;
    vnc_idx = vnc_start_idx:vnc_end_idx ; 
else
    vnc_seg = [] ; 
    vnc_idx = [] ; 
end
end

% -------------------------------------------------------------------------
%% look for neuron number 
% (e.g. if there are multiple segmented neurons from the same hemilineage)
function neuron_num = get_neuron_num(hl_str, delim)
neuron_num_idx = regexp(hl_str, [delim '\d' delim]) ;
if isempty(neuron_num_idx)
    neuron_num = [] ; 
elseif numel(neuron_num_idx) == 1
    neuron_num = hl_str(neuron_num_idx+1) ; 
else
    disp('unsure of what to do about neuron number')
    keyboard
end
end

% -------------------------------------------------------------------------
%% look for driver name 
%(in the form SSXXXXX where X is an integer)
function driver_name = get_driver_name(hl_str)
[driver_start_idx, driver_end_idx] = regexp(hl_str, 'SS\d\d\d\d\d') ;
if isempty(driver_start_idx)
    driver_name = [] ; 
else
    driver_name = hl_str(driver_start_idx:driver_end_idx) ;
end
end

% -------------------------------------------------------------------------
%% look for fly sex (should be M or F)
function sex = get_fly_sex(hl_str)
% male_idx = regexp(hl_str, '[M]') ;
% female_idx = regexp(hl_str, '[F]') ;
if ~contains(hl_str,'male') && ~contains(hl_str,'female')
    sex = [] ; 
elseif contains(hl_str,'male') && ~contains(hl_str,'female')
    sex = 'M' ; 
elseif ~contains(hl_str,'male') && contains(hl_str,'female')
    sex = 'F' ; 
else
    sex = [] ; 
end
end