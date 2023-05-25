% -------------------------------------------------------------------------
% function to load neuron data -- really annoying to have to constantly
% rewrite code to do this
% -------------------------------------------------------------------------
function [neuronData, neuronType] = myLoadNeuronData(neuron_name, ...
    dataType, neuronType, rootPath)
% -------------------
%% inputs and params
if ~exist('dataType','var') || isempty(dataType)
    dataType = 'coords_sym' ;  % ('bw' | 'coords' | 'skel' | 'mask_bw' | 'mask_coords' ) + ('sym' | 'non_sym')
end
if ~exist('rootPath','var') || isempty(rootPath)
    [mfilePath, ~, ~] = fileparts(mfilename('fullpath')) ; 
    parentDirectory = fileparts(mfilePath) ; 
    rootPath = fullfile(parentDirectory, 'data', 'processed_anatomy_data') ;
end
if ~exist('neuronType','var') || isempty(neuronType)
    neuronType = [] ;
end

% general path params
neuronTypeList = {'IN', 'DN', 'MN', 'VUM'} ;
dataPathList = cellfun(@(y) fullfile(rootPath, y), neuronTypeList, ...
    'UniformOutput', false) ;
fileExtList = {'.nrrd', '.tif'} ; % going to search based on raw image


% ---------------------------------------------------------------------
%% search through directories for neuron with given name
% (if we're not given this as input)
if ~isempty(neuronType)
    dataPathMatch = fullfile(rootPath, neuronType) ;
else
    dataPathMatch = [] ;
    % loop over neuron type directories
    for k = 1:length(dataPathList)
        % current neuron directory
        dataPathCurr = dataPathList{k} ;
        
        % for loop over image file formats
        for n = 1:length(fileExtList)
            searchPathCurr = fullfile(dataPathCurr, ['*' fileExtList{n}]) ;
            dirCurr = dir(searchPathCurr) ;
            
            N_files = length(dirCurr) ;
            % names_curr = cell(N_files,1) ;
            
            % try to find file corresponding to neuron in this directory
            for m = 1:N_files
                [~, fn_no_ext, ~] = fileparts(dirCurr(m).name) ; 
                fn_split = strsplit(fn_no_ext,'_') ;
                
                % if we find out neuron, assign the data path
                if strcmpi(fn_split{1},neuron_name) 
                    dataPathMatch = dataPathCurr ;
                    neuronType = neuronTypeList{k} ; 
                    break ;
                end
            end
        end
    end
end

% ---------------------------------------------------------------
%% search within dataPathMatch (if we foudn one)
% first check that we found a matching directory
if isempty(dataPathMatch)
    fprintf('Error finding neuron %s \n', neuron_name)
    keyboard
end

% if we did find a match, look in the appropriate folder (determined by
% dataType)
switch dataType
    case 'coords_sym'
        dataFolder = fullfile(dataPathMatch, 'binarized_sym') ; 
        searchSuffix = '*_coords.mat' ; 
    case 'coords_non_sym'
        dataFolder = fullfile(dataPathMatch, 'binarized_non_sym') ; 
        searchSuffix = '*_coords.mat' ; 
    case 'skel_sym' 
        dataFolder = fullfile(dataPathMatch, 'skeletonized_sym') ; 
        searchSuffix = '*_skel_struct.mat' ; 
    case 'skel_non_sym'
        dataFolder = fullfile(dataPathMatch, 'skeletonized_non_sym') ;
        searchSuffix =  '*_skel_struct.mat' ;
    case {'mask_bw_sym', 'mask_bw_non_sym'} 
        dataFolder = fullfile(dataPathMatch, 'io_masks') ;
        searchSuffix = '*_mask.mat' ;
    case {'mask_coords_sym', 'mask_coords_non_sym'}
        dataFolder = fullfile(dataPathMatch, 'io_masks') ;
        searchSuffix = '*_mask_coords.mat' ;
    otherwise
        fprintf('Error: invalid data type: %s \n', dataType)
end

% get directory of (hopefully) the correct data folder
dataDir = dir(fullfile(dataFolder, searchSuffix)) ;

% look through dir to find match of neuron name
N_files = length(dataDir) ;
names_curr = cell(N_files,1) ; 
for q = 1:N_files
    fn_split = strsplit(dataDir(q).name, '_') ; 
    names_curr{q} = fn_split{1} ; 
end

% --------------------------------------------------------------
%% (hopefully) find match and load data
match_idx = cellfun(@(y) strcmpi(neuron_name, y), names_curr) ; 

if sum(match_idx) ~= 1
   fprintf('Error: could not find unique match for %s \n', neuron_name) 
   neuronData = [] ;
else
    neuronDataPath = fullfile(dataDir(match_idx).folder, ...
        dataDir(match_idx).name) ;
    neuronData = importdata(neuronDataPath) ; 
end

end