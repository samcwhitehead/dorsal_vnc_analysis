% -------------------------------------------------------------------------
% Script to loop over data files (.mat files resulting from
% kinefly_analysis.m) in a given directory and generate combined data
% structures for these multiple files. It is assumed that data are grouped
% by genotype or some other experimental condition
% -------------------------------------------------------------------------
%% PARAMS AND PATH INFO
% Provide the path to look for files within (dataRoot). this code assumes
% that data is further grouped within the supplied folder by experiment 
% type/condition. This would look something like:
% 
% dataRoot/
%   - closed_loop_50x/
%       - 20170828T150500_SS01062
%           - 20170828_SS01062_0000_axo_processed.mat
%       - 20170828T151600_SS01062
%           - ...
%       - ...
%   - closed_loop_20x/
%       - ...
%   - ...
% 
[mfilePath, ~, ~] = fileparts(mfilename('fullpath')) ; 
parentDirectory = fileparts(mfilePath) ;
dataRoot = fullfile(parentDirectory, 'data') ;
savePath = fullfile(parentDirectory,'figs','output') ;

% find all data files within this directory
dataDir = dir(fullfile(dataRoot, '**', '**', '*axo_processed.mat')) ; 

% determine the unique data folders within this directory
dataPaths = arrayfun(@(x) fileparts(x.folder), dataDir, 'UniformOutput',0);
dataPaths = unique(dataPaths) ; 

% save resulting data structures?
saveFlag = true ; 

% ------------------------------------------------------------
%% AGGREGATE DATA INTO STRUCTURES AT FULL RESOLUTION
% loop over unique dataPaths and create 'time_series_struct' for the data
% contained in each folder
for k = 1:length(dataPaths)
    % current data path
    dataPath = dataPaths{k} ;
    
    % pull out 'condition' from this folder name
    [~, condition, ~] = fileparts(dataPath) ; 
    
    % get time series structure
    time_series_struct = get_time_series_struct(dataPath, condition,...
        saveFlag);
    
    fprintf('Completed time series struct for: \n %s \n', dataPath) 
end

% ------------------------------------------------------------
%% CREATE GRAND MEAN STRUCT ACROSS CONDITIONS

% generate grand_mean_struct
grand_mean_struct = laptop_get_grand_mean_struct(dataPaths) ;

% save output?
if saveFlag
    save(fullfile(dataRoot, 'grand_mean_struct.mat'),'grand_mean_struct')
end

fprintf('Completed grand mean struct for: \n %s \n', dataRoot) 