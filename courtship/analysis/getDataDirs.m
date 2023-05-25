%--------------------------------------------------------------------------
% much of the analysis infrasructure is already in place from code base.
% this script just grabs the paths/genotypes for segmented, analyzed song
% data, so that we can run full analysis
%--------------------------------------------------------------------------
%% PATH INFO
% path to data file 
[mfilePath, ~, ~] = fileparts(mfilename('fullpath')) ; 
parentDirectory = fileparts(mfilePath) ;
rootPath = fullfile(parentDirectory, 'data', 'example') ;
savePath = rootPath ;

dataPathStruct = struct() ;

saveFlag = true ; 
LLR_threshold = 0 ; 

% all folders that have already been run through BatchFlySongAnalsis.m
dataDir = dir(fullfile(rootPath, ...
    ['*Results_LLR=' num2str(LLR_threshold) '_*'])) ; 
dataDir = dataDir([dataDir.isdir]) ; 

% PlotAnalysis takes in a cell array of folders
StatsPathCell = arrayfun(@(x) fullfile(x.folder,x.name),dataDir,...
    'UniformOutput',false) ; 

% get folder locations for un-batch-analyzed data
DataFolderSplit = arrayfun(@(x) strsplit(x.name,'_'), dataDir,...
    'UniformOutput',false) ; 
DataFolderPrefixes = cellfun(@(x) x(1), DataFolderSplit) ; 
DataPathCell = arrayfun(@(x,y) fullfile(x.folder,strcat(y,'_out')), dataDir,...
    DataFolderPrefixes) ; %, 'UniformOutput',false) ;

validControlTypes = {'SS01062','SS01055'} ;

%--------------------------------------------------------------------------
%% pull out genotypes contained in each 
GenotypeCell = cell(length(StatsPathCell),1) ; 
GenotypeChannelCell = cell(length(StatsPathCell),1) ; 
GenotypeTypeCell = cell(length(StatsPathCell),1) ; % couldn't resist...
GenotypeTypeCell_valid = cell(length(StatsPathCell),1) ; 
for i = 1:length(GenotypeCell)
    folderDir = dir([StatsPathCell{i} '\PS*']) ; 
    
    % get datetime at which analysis occured
    folderParentPath = folderDir(1).folder ; 
    [~,folderParent,~] = fileparts(folderParentPath) ; 
    folderParent_split = strsplit(folderParent,'_') ;
    analysisDate = folderParent_split{end} ; 
    
    genotype_curr = cell(1,length(folderDir)) ; 
    channel_num_list = zeros(1,length(folderDir)) ; 
    for j = 1:length(folderDir) 
        fileName = folderDir(j).name ; 
        fileName_split = strsplit(fileName,'_') ; 
        
        % it looks like genotye is in filename between channel number and
        % datetime for analysis--going to try to grab it out
        channel_ind = find(contains(fileName_split,'ch')) ; 
        analysis_date_ind = find(contains(fileName_split,analysisDate)) ; 
        
        genotype_curr{j} = ...
            strjoin(fileName_split((channel_ind+1):(analysis_date_ind - 1)),'_') ;
        
        channel_num = str2double(fileName_split{channel_ind}(3:end)) ;
        channel_num_list(j) = channel_num ; 
    end
    
    [genotype_unique, ~, geno_channel_idx] = unique(genotype_curr) ; 
    channel_for_genotype = cell(1,length(genotype_unique)) ; 
    for k = 1:length(genotype_unique)
       channel_for_genotype{k} = channel_num_list(geno_channel_idx == k) ; 
    end
    
    GenotypeCell{i} = genotype_unique ;
    GenotypeTypeCell{i} = cellfun(@(x) contains(x, 'control'),genotype_unique) ; 
    GenotypeTypeCell_valid{i} = cellfun(@(x) contains(x, 'control') & ...
        any(contains(x,validControlTypes)),genotype_unique) ; 
    GenotypeChannelCell{i} = channel_for_genotype ; 
end

% in case we can run analysis in one function call, get all unique
% genotypes and whether or not they are control

GenotypeCellAll = [GenotypeCell{:}] ; 
GenotypeTypeCellAll = [GenotypeTypeCell{:}] ; 
GenotypeTypeCellAll_valid = [GenotypeTypeCell_valid{:}] ; 

[GenotypeCellAll_unique, uniqueInd,~]  = unique(GenotypeCellAll) ; 

%uniqueInd_logical = false(size(GenotypeCellAll)) ; 
%uniqueInd_logical(uniqueInd) = true ; 

%validControlInd = cellfun(@(x) any(contains(x,validControlTypes)),GenotypeCellAll) ; 
GenotypeTypeCellAll_unique = GenotypeTypeCellAll(uniqueInd) ; 
GenotypeTypeCellAll_unique_valid = GenotypeTypeCellAll_valid(uniqueInd) ; 


%--------------------------------------------------------------------------
%% add fields to struct
dataPathStruct.rootPath = rootPath ; 
dataPathStruct.FoldersStats = StatsPathCell ;
dataPathStruct.FoldersData = DataPathCell ; 
dataPathStruct.GenotypeByFolder = GenotypeCell ;
dataPathStruct.ControlIndByFolder =  GenotypeTypeCell ; 
dataPathStruct.ControlIndByFolder_valid =  GenotypeTypeCell_valid; 
dataPathStruct.GenotypeChannelByFolder = GenotypeChannelCell ;

dataPathStruct.GenotypeAll = GenotypeCellAll_unique ; 
dataPathStruct.ControlIndAll = GenotypeTypeCellAll_unique ; % true if control, false if experimental
dataPathStruct.ControlIndAll_valid = GenotypeTypeCellAll_unique_valid ; % true if control and one of the blank split lines, false if experimental

%--------------------------------------------------------------------------
%% save results?
if saveFlag 
   save(fullfile(savePath, 'dataPathStruct.mat'),'dataPathStruct') 
    
end
