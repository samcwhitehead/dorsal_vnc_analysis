% -------------------------------------------------------------------------
% Script to generate "songDataStruct", which is used to make main song
% stats plot and contains aggregated pulse/sine song data across genotypes
% (obtained from getSongStats.m)
%
% -------------------------------------------------------------------------
% load structure with data directory information
[mfilePath, ~, ~] = fileparts(mfilename('fullpath')) ; 
parentDirectory = fileparts(mfilePath) ;
rootPath = fullfile(parentDirectory, 'data', 'example') ;
savePath = rootPath ;
dataDirName = 'dataPathStruct.mat' ;
dataPathStruct = importdata(fullfile(rootPath, dataDirName)) ;

% get all folders for genotypes
dataFolder = fullfile(rootPath, 'frac_singing') ;
dataDir = dir(fullfile(dataFolder, '*')) ;
dataDir = dataDir(3:end) ;

% define the lines/variables that we want to plot (currently limited by
% what's in data structure
MN_cell = {'DVM', 'tp2' ,  'ps1', 'i2' , 'hg1', 'hg2', 'hg3'} ; % , 'hg1 + hg3/4', 'b1 + hg4?', 'hg1 + i1', 'i1 + i2'} ; % which motor neuron drivers to plot
driver_cell = {'SS41068', 'SS47120', 'SS47152', 'SS37246', 'SS48311',...
    'SS37253','SS49039'} ; % , 'SS40864', 'SS40980', 'SS45772', 'SS45782'} ; % for the moment, pick out good examples of each MN
var_cell = {'fracSinging','pulseIndex','sineIndex'} ;
varName_cell = {'% Singing', 'Pulse Index', 'Sine Index'} ;

silencer_cell = {'Kir'} ;
exprCond_cell = {'ctrl','expr'} ;

define_constants_kf ;

% save results?
if ~exist(savePath,'dir')
    mkdir(savePath)
end
saveFlag = false ;

% params
singFracMin = 5e-3 ; %5e-3 ;
singPercThresh = singFracMin ; % 5e-3 ;
N_plot_rows = length(var_cell) ;
N_plot_cols = length(driver_cell) ;


%--------------------------------------------------------------------------
%% for convenience/reproducibility, gather data to be plotted into struct
GenotypeAll = dataPathStruct.GenotypeAll ;
ControlIndAll = dataPathStruct.ControlIndAll;
GenotypeExpr = GenotypeAll(~ControlIndAll) ;

songDataStruct = struct() ;
%cc = 1 ;

for i = 1:length(driver_cell)
    driver_ind = arrayfun(@(x) contains(x.name,driver_cell{i}),dataDir) ;
   
    for j = 1:length(silencer_cell)
        silencer_ind =  arrayfun(@(x) contains(x.name,silencer_cell{j}),...
            dataDir) ;
        curr_folder_ind = driver_ind & silencer_ind ;
        
        if sum(curr_folder_ind) ~= 1
            disp('could not find data for this genotype -- skipping')
            continue
        end
        
        dataPath_curr = fullfile(dataDir(curr_folder_ind).folder,...
            dataDir(curr_folder_ind).name,'songFracStruct.mat') ;
        songFracStruct = importdata(dataPath_curr) ;
        
        % store some general info about the fly line
        songDataStruct(i).driver = driver_cell{i} ;
        songDataStruct(i).silencerType = silencer_cell{j} ;
        songDataStruct(i).dataPath = dataPath_curr;
        songDataStruct(i).ExprGenotype = GenotypeExpr{curr_folder_ind};
        %--------------------------------------------------------------------------
        %% identify flies that sing for less that a given fraction of the recording
        song_frac = [songFracStruct(:).songIndex] ;
        good_ind = (song_frac >= singFracMin) ;
        singing_ind = (song_frac > singPercThresh) ;
        
        ctrl_ind = arrayfun(@(x) strcmp(x.flyType,'ctrl'),songFracStruct) ;
        expr_ind = arrayfun(@(x) strcmp(x.flyType,'expr'),songFracStruct) ;
        
        songDataStruct(i).([var_cell{1} '_' silencer_cell{j} '_ctrl']) = ...
            100*sum(singing_ind & ctrl_ind)/ sum(ctrl_ind) ;
        songDataStruct(i).([var_cell{1} '_' silencer_cell{j} '_expr']) = ...
            100*sum(singing_ind & expr_ind)/ sum(expr_ind) ;
        
        songDataStruct(i).SongFracTable = table([sum(singing_ind & expr_ind) ; ...
            sum(~singing_ind & expr_ind)], [sum(singing_ind & ctrl_ind) ; ...
            sum(~singing_ind & ctrl_ind)],'VariableNames',{'Expr','Ctrl'},...
            'RowNames',{'Singing','NoSinging'}) ;
        
        try
            songFracStruct = songFracStruct(good_ind) ;
        catch
            disp('no good trials')
            for k = 2:length(var_cell)
                songDataStruct(i).([var_cell{k} '_' silencer_cell{j} ...
                    '_ctrl']) = [] ;
                songDataStruct(i).([var_cell{k} '_' silencer_cell{j} ...
                    '_expr']) = [] ;
            end
            %cc = cc + 1 ;
            continue
        end
        
        ctrl_ind = arrayfun(@(x) strcmp(x.flyType,'ctrl'),songFracStruct) ;
        expr_ind = arrayfun(@(x) strcmp(x.flyType,'expr'),songFracStruct) ;
        ctrl_sub = find(ctrl_ind) ;
        expr_sub = find(expr_ind) ;
        
        for k = 2:length(var_cell)
            songDataStruct(i).([var_cell{k} '_' silencer_cell{j} ...
                '_ctrl']) = [songFracStruct(ctrl_ind).(var_cell{k})] ;
            songDataStruct(i).([var_cell{k} '_' silencer_cell{j} ...
                '_expr']) = [songFracStruct(expr_ind).(var_cell{k})] ;
        end
        
    end
end

% ---------------------------------------------
%% save results?
if saveFlag
   save(fullfile(savePath, 'songDataStruct.mat'), 'songDataStruct') 
end