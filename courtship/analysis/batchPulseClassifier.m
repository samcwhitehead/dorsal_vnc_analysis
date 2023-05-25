%--------------------------------------------------------------------------
% a script to run murthy lab pulse classification across all genotypes
%--------------------------------------------------------------------------
%% PATH AND PARAMS
% load structure with data directory information
% load structure with data directory information
[mfilePath, ~, ~] = fileparts(mfilename('fullpath')) ; 
parentDirectory = fileparts(mfilePath) ;
% path to data directory struct created by getDataDirs.m:
rootPath = fullfile(parentDirectory, 'data', 'example') ;
% path to folder for saving pulseClassStruct (output)
savePath = fullfile(rootPath, 'pulse_classifier') ;

% load data directory structure
dataDirName = 'dataPathStruct.mat' ;
dataPathStruct = importdata(fullfile(rootPath, dataDirName)) ;

% save results?
saveFlag = true ;

% plot classifier results along the way?
plotFlag = true ; 

% define params
params = loadSongParams() ; 

%% get info on all genotypes so they can be looped through
GenotypeAll = dataPathStruct.GenotypeAll ;
ControlIndAll = dataPathStruct.ControlIndAll;
ControlIndAll_valid = dataPathStruct.ControlIndAll_valid ;

dataFolders = dataPathStruct.FoldersData ;
GenotypeChannels = dataPathStruct.GenotypeChannelByFolder ;
GenotypeByFolder = dataPathStruct.GenotypeByFolder ;
ControlIndByFolder_valid = dataPathStruct.ControlIndByFolder_valid ;

% get full list of experimental genotypes
GenotypeExpr = GenotypeAll(~ControlIndAll) ;

%% loop through genotypes to pull data from all experiments/analyze
for i = 1:length(GenotypeExpr)
    tic
    ExprGenotype = GenotypeExpr{i} ;
    dataFolderInd = find(cellfun(@(x) any(strcmp(x,ExprGenotype)),...
        GenotypeByFolder)) ;
    
    pulseClassStruct = struct() ;
    cc = 1 ;
    if saveFlag 
        savePath_curr = fullfile(savePath, ExprGenotype) ;
        if ~exist(savePath_curr,'dir')
            mkdir(savePath_curr)
        end
    end
    % loop through all folders containing given experimental genotype
    for j = 1:length(dataFolderInd)
        dataFolderCurr = dataFolders{dataFolderInd(j)} ;
        GenotypesInFolder = [GenotypeByFolder{dataFolderInd(j)}] ;
        GenotypeChannelsCurr = [GenotypeChannels{dataFolderInd(j)}] ;
        
        % find data files corresponding to correct genotypes (valid
        % controls, etc.)
        ControlIndInFolder = find([ControlIndByFolder_valid{dataFolderInd(j)}]) ;
        ExprIndInFolder = find(strcmp(GenotypesInFolder,ExprGenotype)) ;
        GenotypesCurrInd = [ExprIndInFolder, ControlIndInFolder] ;
        %GenotypesCurr = GenotypesInFolder(GenotypesCurrInd) ;
        
        % analyze each channel
        for k = GenotypesCurrInd
            ChannelsCurr = GenotypeChannelsCurr{k} ;
            for m = 1:length(ChannelsCurr)
                channel_num = ChannelsCurr(m) ;
                dataDir = rdir(fullfile(dataFolderCurr, ['\*ch' ...
                    num2str(channel_num) '.mat'])) ;
                if length(dataDir) ~= 1
                    disp('error loading data')
                    continue
                else
                    dataPath = dataDir(1).name ;
                end
                
                % load data analyzed by BatchSongAnalysis.m
                data_in = load(dataPath) ;
                
                % perform pulse classification using Murthy Lab classifier
                try
                    [pulsesToAnalyze, pulseMat, pulsesNorm, ...
                        pulseLabels, h_pclass] = ...
                        myPulseClassifier(data_in, params, plotFlag) ;
                catch
                    disp('Error determining pulse classification')
                    continue
                end
                
                % add fields to structure
                pulseClassStruct(cc).genotype = GenotypesInFolder{k} ;
                pulseClassStruct(cc).dataPath = dataPath ;
                pulseClassStruct(cc).pulsesToAnalyze = pulsesToAnalyze ;
                pulseClassStruct(cc).pulseMat = pulseMat ;
                pulseClassStruct(cc).pulsesNorm = pulsesNorm ;
                pulseClassStruct(cc).pulseLabels = pulseLabels ;
                if strcmp(GenotypesInFolder{k}, ExprGenotype)
                    pulseClassStruct(cc).flyType = 'expr' ; 
                else
                    pulseClassStruct(cc).flyType = 'ctrl' ;
                end
                pulseClassStruct(cc).params = params ;
                
                cc = cc + 1 ; 
                
                if saveFlag && plotFlag
                   savefig(h_pclass,fullfile(savePath_curr,...
                       [GenotypesInFolder{k} '_ch' num2str(channel_num) '.fig'])) 
                   print(h_pclass,fullfile(savePath_curr,...
                       [GenotypesInFolder{k} '_ch' num2str(channel_num) '.png']),...
                       '-dpng','-r300')
                   
                   close all
                end
            end
            
        end
    end
    
    if saveFlag
        save(fullfile(savePath_curr, 'pulseClassStruct.mat'),...
            'pulseClassStruct')
    end
    toc
    disp(['Completed ' num2str(i) '/' num2str(length(GenotypeExpr))])
end