%--------------------------------------------------------------------------
% function to generate time series structures for a set of kinefly data
% folders. these structures just compile all of the data sets in one folder
%--------------------------------------------------------------------------
function time_series_struct = ...
    get_time_series_struct(dataPath,condition, saveFlag)
%--------------------------------------------------------------------------
%% params
define_constants_kf

% which variables to grab from the axo_processed folders. variables like
% 'R_AMP' are defined in 'define_constants_kf.m'
var_oi = [R_AMP, L_AMP, R_DEV_1, L_DEV_2,...
    R_DEV_2, L_DEV_1, WBF] ;
var_str = {'R_AMP', 'L_AMP', 'R_DEV_1', 'L_DEV_2',...
    'R_DEV_2', 'L_DEV_1', 'WBF'} ;

%% get data directory
% within dataPath, there should be folders that each contain a single
% kinefly run. we want to first get the directory for all of these
dataDir = dir([dataPath '\20*']) ;

%% generate structure
time_series_struct = struct('driver',[]) ;
%--------------------------------------------------------------------------
% pull data from individual files
for i = 1:length(dataDir)
    %get fly genotype
    folderName = dataDir(i).name ;
    folderName_split = strsplit(folderName,'_') ;
    driver = folderName_split{2} ;
    time_series_struct(i).driver = driver ;
    folderPath = fullfile(dataDir(i).folder, dataDir(i).name) ; 
    
    %get time series data. subsample?
    if strcmp(condition,'open_loop')
        fileDir = dir(fullfile(folderPath, ['\*_' open_loop_type ...
            '_axo_processed.mat'])) ;
    else
        fileDir = dir(fullfile(folderPath, '\*axo_processed.mat')) ;
    end
    
    try
        data = importdata(fullfile(fileDir(1).folder, fileDir(1).name)) ;
    catch
        continue
    end
    
    for j = 1:length(var_str)
        var_curr = var_str{j} ; 
        time_series_struct(i).(var_curr) = data.output{var_oi(j)} ;
    end
end

%% remove any empty entries
empty_ind = arrayfun(@(x) any( structfun(@isempty, x) ),time_series_struct) ;
time_series_struct = time_series_struct(~empty_ind) ; 

%% save results?
if saveFlag
    if strcmp(condition,'open_loop')
        save([dataPath '\time_series_struct_' open_loop_type '.mat'],...
            'time_series_struct','-v7.3')
    else
        save([dataPath '\time_series_struct.mat'],'time_series_struct','-v7.3')
    end
end


end