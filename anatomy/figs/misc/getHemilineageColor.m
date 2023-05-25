% -------------------------------------------------------------------------
% function to keep track of hemilineage color scheme
%
% INPUTS:
%   - hemilineageIDs: cell array of strings (or one string) that contain
%       neuroblast number and notch type, e.g. {'18B', '0A'} | '18B'
%   - colorMode: string defining how we should color hemilineages (each
%       hemilineage gets its own color, or by behavior/neurotransmitter)
% -------------------------------------------------------------------------
function [colorsOut, colorKey] = getHemilineageColor(hemilineageID,...
    colorMode)
% ----------------------------
% inputs
if ~exist('colorMode','var') || isempty(colorMode)
    colorMode = 'hemilineage' ; % 'hemilineage' | 'behavior' | 'neurotransmitter'
end
% ---------------------------------------------
% list hemilineage types in our collection
hemilineage_names = {'0A', '0B', '2A', '3B', '5B', '6A', '6B', '7B', ...
    '8B', '11A', '11B', '12A', '17A', '18B', '19A', '19B'} ;

% also list names for embryonic and abdominal lineages
other_lin_names = {'abd', 'emb'} ;

% combine to one cell
names_all = [hemilineage_names, other_lin_names] ;

% make sure hemilineageID is a cell so we can use cellfun later on
if ~iscell(hemilineageID)
    hemilineageID = {hemilineage_ID} ;
end
% ------------------------------------------------
% get colormap for all hemilineages
N_hl = length(hemilineage_names) ;

% load additional data if required
if strcmpi(colorMode, 'behavior') || strcmpi(colorMode, 'neurotransmitter')
    hemilineage_table = load_hemilineage_info() ;
end

% colormap depends on how what we're using as basis for colors
switch colorMode
    case 'hemilineage'
        colorsHL = viridis(N_hl) ;
        colorsHL = flipud(colorsHL) ; % make earlier hemilineages light and latter dark
        
        % also add colors for other lineage types -- use gray
        colorsOther = gray(length(other_lin_names) + 2) ;
        colorsOther = colorsOther(2:end-1,:) ; % remove pure white and black
        
        % combine color matrices
        colorsAll = [colorsHL ; colorsOther ] ;
        
        % list name corresponding to each color
        colorKey = names_all' ;
        
    case 'behavior'
        % get all behaviors listed in hemilineage table
        behaviors_all = hemilineage_table.AssociatedBehaviors ;
        
        % since some rows contain multiple entries, need to break them up
        behaviors_all = cellfun(@(y) strsplit(y,', '), behaviors_all, ...
            'UniformOutput', false) ;
        behaviors_all = horzcat(behaviors_all{:}) ;
        
        % get unique behaviors from this
        behaviors_unique = unique(behaviors_all) ;
        
        % create color matrix with an entry for each behavior (including
        % for blank entries ("-"), which will be gray)
        colorsBehavior = vertcat(0.7*[1,1,1], ...
            brewermap(length(behaviors_unique) - 1, 'Set2')) ;
        
        % read out main behaviors from table
        behaviors_main = hemilineage_table.MainBehavior ;
        
        % provide names for each color
        colorKeyTemp = behaviors_unique ;
        % ... but replace "-" with unknown
        colorKeyTemp{1} = 'unknown' ;
        colorKey = cell(length(behaviors_main),1) ;
        
        % now use color matrix above to make another -- this one with a row
        % for each hemilineage. we can index on this later
        colorsAll = zeros(length(behaviors_main), 3) ;
        
        % loop over hemilineage table rows
        for q = 1:length(behaviors_main)
            % find index in list of unique behavior that matches behavior
            % of current row (row = hemilineage)
            match_idx = cellfun(@(y) strcmpi(y, behaviors_main{q}), ...
                behaviors_unique) ;
            
            % make sure we find a match
            if sum(match_idx) ~= 1
                fprintf('Error finding match for behavior %s \n', ...
                    behaviors_main{q})
                keyboard
            end
            
            % assign color and key entry for this row
            colorsAll(q,:) = colorsBehavior(match_idx,:) ;
            colorKey{q} = colorKeyTemp{match_idx} ;
        end
        
        
        
    case 'neurotransmitter'
        % get all neurotransmitters listed in hemilineage table
        neurotrans_all = hemilineage_table.Neurotransmitter ;
        
        % get unique neurotransmitters
        [neurotrans_unique, ~, ic] = unique(neurotrans_all) ;
        
        % create color matrix with an entry for each neurotransmitter
        % (including for blank entries ("-"), which will be gray)
        colorsSet3 = brewermap(length(neurotrans_unique) + 2, 'Set3') ; % tweak to avoid pastel yellow
        colorsNeurotrans = vertcat(0.7*[1,1,1], colorsSet3(4:end,:)) ;
        
        % now use color matrix above to make another -- this one with a row
        % for each hemilineage. we can index on this later
        colorsAll = colorsNeurotrans(ic,:) ;
        
        % -------------------------------------
        % give names for colors in this case
        colorKey = neurotrans_unique ;
        colorKey{1} = 'unknown' ; % replace "-" with "unknown"
        
        % make it so that we have a key entry for each row in "colorsAll"
        colorKey = colorKey(ic) ;
    otherwise
        fprintf('Invalid colorMode selection: %s \n', colorMode)
        keyboard
end

% -------------------------------------------------
% find indices in "colorsAll" for inputs
color_ind = cellfun(@(y) find(strcmpi(y,names_all),1, 'first'), ...
    hemilineageID) ;


% ---------------------------------------------------
% grab just colors/keys indicated by index
colorsOut = colorsAll(color_ind,:) ;
colorKey = colorKey(color_ind) ;

end