%--------------------------------------------------------------------------
% a script to grab fraction of singers, sine index, and pulse index from a
% group of flies
%--------------------------------------------------------------------------

% load structure with data directory information
[mfilePath, ~, ~] = fileparts(mfilename('fullpath')) ; 
parentDirectory = fileparts(mfilePath) ;
rootPath = fullfile(parentDirectory, 'data', 'example') ;
dataDirName = 'dataPathStruct.mat' ;
dataPathStruct = importdata(fullfile(rootPath, dataDirName)) ;

savePath = fullfile(rootPath, 'frac_singing') ;

saveFlag = true ;
plotFlag = false ;

Nbins = 25 ;
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
    
    songFracStruct = struct() ;
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
                dataDir = rdir(fullfile(dataFolderCurr, ...
                    ['\*ch' num2str(channel_num) '.mat'])) ;
                if length(dataDir) ~= 1
                    disp('error loading data')
                    continue
                else
                    dataPath = dataDir(1).name ;
                end
                
                data_in = load(dataPath) ;
                
                try
                    songIndexStruct = ...
                        getPulseAndSineIndices(data_in, params) ;
                catch
                    disp('failed to get pulse/sine indices')
                    continue
                end
                
                songFracStruct(cc).genotype = GenotypesInFolder{k} ;
                songFracStruct(cc).dataPath = dataPath ;
                
                structFieldnames = fieldnames(songIndexStruct) ; 
                for fn_ind = 1:length(structFieldnames)
                    songFracStruct(cc).(structFieldnames{fn_ind}) = ...
                        songIndexStruct.(structFieldnames{fn_ind}) ; 
                end
                songFracStruct(cc).Fs = data_in.Data.fs ;
                songFracStruct(cc).recording_duration = ...
                    length(data_in.Data.d) ;
                if strcmp(GenotypesInFolder{k}, ExprGenotype)
                    songFracStruct(cc).flyType = 'expr' ;
                else
                    songFracStruct(cc).flyType = 'ctrl' ;
                end
                songFracStruct(cc).params = params ;
                
                cc = cc + 1 ;
                
            end
            
        end
    end
    
    if plotFlag
        h_hist = figure('PaperPositionMode','auto') ;
        ctrl_ind = arrayfun(@(x) strcmp(x.flyType, 'ctrl'), songFracStruct) ;
        histogram([songFracStruct(ctrl_ind).songIndex],Nbins,...
            'normalization','count','DisplayStyle','stairs','EdgeColor','b')
        hold on
        histogram([songFracStruct(~ctrl_ind).songIndex],Nbins,...
            'normalization','count','DisplayStyle','stairs','EdgeColor','r')
        title(ExprGenotype)
        legend({'Ctrl','Expr'})
    end
    if saveFlag
        save(fullfile(savePath_curr, 'songFracStruct.mat'),...
            'songFracStruct')
        if plotFlag
            print(h_hist,fullfile(savePath_curr,...
                [GenotypesInFolder{k} 'frac_singing.png']),...
                '-dpng','-r300')
            
            close all
        end
    end
    toc
    disp(['Completed ' num2str(i) '/' num2str(length(GenotypeExpr))])
end