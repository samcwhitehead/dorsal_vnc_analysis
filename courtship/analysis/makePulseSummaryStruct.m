% -------------------------------------------------------------------------
% Script to generate a combined data structure called "pulseSummaryStruct"
% that incorporates data across multiple "pulseClassStruct" (the output of
% batchPulseClassifier.m that contains data for a single genotype).
%
% This combined data structure is used to generate the pulse classification
% figure (see plotPulseFig.m)
%
% -------------------------------------------------------------------------
%% PATH AND PARAMS
rootPath = '' ;   % path to 'dataPathStruct.mat'
savePath = '' ;   % path to folder where summary struct is saved
pulsePath = '' ;  % path to fldr containing "pulseClassStruct" per genotype
dataDirName = 'dataPathStruct.mat' ;
dataPathStruct = importdata(fullfile(rootPath, dataDirName)) ;

% define the lines/variables that we want to plot (currently limited by
% what's in data structure
MN_cell = {'DVM', 'tp2' ,  'ps1', 'i2' , 'hg1', 'hg2', 'hg3'} ;  % which motor neuron drivers to plot
driver_cell = {'SS41068', 'SS47120', 'SS47152', 'SS37246', 'SS48311',...
    'SS37253','SS49039'} ; % for the moment, pick out good examples of each MN
var_cell = {'pulseFrac','slow','fast'} ;
varName_cell = {'P_{slow} frac.', 'Slow\nPulse', 'Fast\nPulse'} ;

silencer_cell = {'Kir'} ;
exprCond_cell = {'ctrl','expr'} ;

% save result?
if ~exist(savePath,'dir')
    mkdir(savePath)
end
saveFlag = false ;

% -------------------------------
% general params
define_constants_kf ;

singFracMin = 5e-3 ;
singPercThresh = 5e-3 ;
minPulseNum = 20 ;

N_BOOT = 500 ;
Fs = 10000 ;
xlim_t = 10*[-1, 1] ;

N_plot_rows = length(var_cell) ;
N_plot_cols = length(driver_cell) ;

%--------------------------------------------------------------------------
%% loop over genotypes, gather all data to be plotted into struct
GenotypeAll = dataPathStruct.GenotypeAll ;
ControlIndAll = dataPathStruct.ControlIndAll;
GenotypeExpr = GenotypeAll(~ControlIndAll) ;

% initialize struct and path to save to
summaryStructPath = fullfile(savePath, 'pulseSummaryStruct.m') ;
pulseSummaryStruct = struct() ;

for i = 1:length(driver_cell)
    for j = 1:length(silencer_cell)
        % get folder containing the correct driver line and silencer type
        driver_ind = cellfun(@(y) contains(y, driver_cell{i}),...
            GenotypeExpr) ;
        silencer_ind = cellfun(@(y) contains(y, silencer_cell{j}),...
            GenotypeExpr) ;
        comb_ind = driver_ind & silencer_ind ;
        
        if sum(comb_ind) ~= 1
            disp('failed to properly specify folder')
            keyboard
        end
        exprGenotype = GenotypeExpr{comb_ind} ;
        
        % load both pulse classification structure and song frac struct
        pulseClassStruct = importdata(fullfile(pulsePath, ...
            exprGenotype, 'pulseClassStruct.mat')) ;
        
        % store some general info about the fly line
        pulseSummaryStruct(i).driver = driver_cell{i} ;
        pulseSummaryStruct(i).silencerType = silencer_cell{j} ;
        pulseSummaryStruct(i).ExprGenotype = exprGenotype ;
        
        %--------------------------------------------------------------------------
        % identify flies that sing for less that a given fraction of the recording
        numPulses = arrayfun(@(x) length(x.pulseLabels),...
            pulseClassStruct) ;
        good_ind = (numPulses > minPulseNum) ;
        pulseClassStruct = pulseClassStruct(good_ind) ;
        
        %------------------------------------------------
        %% first get fraction of slow vs fast pulse mode
        ctrl_ind = arrayfun(@(x) strcmp(x.flyType,'ctrl'),...
            pulseClassStruct) ;
        expr_ind = arrayfun(@(x) strcmp(x.flyType,'expr'),...
            pulseClassStruct) ;
        pulseSlowFrac = ...
            arrayfun(@(x) sum(x.pulseLabels==1)/length(x.pulseLabels),...
            pulseClassStruct) ;
        pulseSummaryStruct(i).([var_cell{1} '_' silencer_cell{j} '_ctrl']) = ...
            pulseSlowFrac(ctrl_ind) ;
        pulseSummaryStruct(i).([var_cell{1} '_' silencer_cell{j} '_expr']) = ...
            pulseSlowFrac(expr_ind) ;
        
        %------------------------------------------------
        %% next get average of fast and slow
        for k = 2:3
            data_curr = arrayfun(@(x) ...
                nanmean(x.pulsesNorm(x.pulseLabels == (k-2),:)), ...
                pulseClassStruct,'UniformOutput',0) ;
            data_curr = cell2mat(data_curr') ;
            pulseSummaryStruct(i).([var_cell{k} '_' silencer_cell{j} ...
                '_ctrl']) = data_curr(ctrl_ind,:) ;
            pulseSummaryStruct(i).([var_cell{k} '_' silencer_cell{j} ...
                '_expr']) = data_curr(expr_ind,:) ;
        end
        
    end
end
clear pulseClassStruct

% -----------------------------------------
%% save result?
if saveFlag
    save(summaryStructPath,'pulseSummaryStruct') ;
end
