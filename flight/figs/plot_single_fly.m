%--------------------------------------------------------------------------
% function to generate plots for the wing angles for a single fly
%--------------------------------------------------------------------------
function h_main = plot_single_fly(dataPath, driver, flyInd, rawTracesFlag,...
    smoothFlag, excludeTrialsFlag)
%--------------------------------------------------------------------------
%% inputs
if ~exist('dataPath','var') || isempty(dataPath)
    [mfilePath, ~, ~] = fileparts(mfilename('fullpath')) ; 
    parentDirectory = fileparts(mfilePath) ;
    dataPath = fullfile(parentDirectory, 'data') ; 
end
if ~exist('rawTracesFlag','var')
    rawTracesFlag = false ;
end
if ~exist('smoothFlag','var') || isempty(smoothFlag)
    smoothFlag = false ;
end
if ~exist('excludeTrialsFlag','var')
    excludeTrialsFlag = false ;
end

ctrl_driver = 'SS01062' ;
flyInd_ctrl = 1 ;
%{
driver = 'SS37253' ;
flyInd = 4 ;
%}
%--------------------------------------------------------------------------
%% plot params
define_constants_kf ;

figPosition = [0.5, 5.5000, 1.55, 4.42] ;
exp_cmap = brewermap(4,cmap_struct.(driver)) ;
expColor = exp_cmap(end,:) ;
ctrl_cmap = brewermap(4,cmap_struct.SS01062) ;
ctrlColor = ctrl_cmap(end,:) ;

lw_thin = 0.5 ;
lw_thick = 1.5 ;
ci_alpha = 0.5 ;
line_alpha = 0.15 ;

xlim = [-0.25, 0.5] ;
%----------------------------------------
% sort out plot types/labels/axis limits
plotVarNames_cell = {{'R_AMP', 'L_AMP'}, ...
    {'R_DEV_1', 'L_DEV_2'}, ...
    {'R_DEV_2', 'L_DEV_1'}, ...
    {'WBF'}} ;
ylim_cell = {12*[-1, 1], 6*[-1, 1], 6*[-1, 1], 9*[-1, 1]} ;
ylabel_cell = {'\Delta Amp (deg)', '\Delta Dev_{fwd} (deg)', ...
    '\Delta Dev_{back} (deg)', '\Delta WBF (Hz)'} ;

%----------------------------------------
% timing info
global_t = T_START:(1/SAMPLE_RATE):T_END ;
sub_time = global_t(1:N_SUB_SAMPLE:end) - LED_start_time; % subtract 1 so LED starts at 0 seconds
LED_start_time = 0 ;
%--------------------------------------------------------------------------
%% load data
dataDir = dir(fullfile(dataPath, 'SS*')) ;
driver_ind = arrayfun(@(x) strcmp(x.name,driver), dataDir) ;
ctrl_driver_ind = arrayfun(@(x) strcmp(x.name,ctrl_driver), dataDir) ;
grand_mean_struct_exp = importdata(fullfile(dataPath, ...
    dataDir(driver_ind).name, 'grand_mean_struct.mat')) ;
grand_mean_struct_ctrl = importdata(fullfile(dataPath, ...
    dataDir(ctrl_driver_ind).name, 'grand_mean_struct.mat')) ;

%--------------------------------------------------------------------------
%% fetch data for plotting
[data_mat_cell_exp, bad_ind_cell_exp] = ...
    fetch_data_mats(grand_mean_struct_exp, flyInd, plotVarNames_cell, ...
    excludeTrialsFlag) ;
[data_mat_cell_ctrl, bad_ind_cell_ctrl] = ...
    fetch_data_mats(grand_mean_struct_ctrl, flyInd_ctrl, plotVarNames_cell,...
    excludeTrialsFlag) ;

%--------------------------------------------------------------------------
%% combine indices for good data traces
bad_ind_exp = false(size(data_mat_cell_exp{1},1), 1) ;
for p = 1:length(bad_ind_cell_exp)
    bad_ind_exp = (bad_ind_exp | bad_ind_cell_exp{p}) ;
end
bad_ind_ctrl = false(size(data_mat_cell_ctrl{1},1), 1) ;
for p = 1:length(bad_ind_cell_ctrl)
    bad_ind_ctrl = (bad_ind_ctrl | bad_ind_cell_ctrl{p}) ;
end

%--------------------------------------------------------------------------
%% generate figure
h_main = figure('PaperPositionMode','auto','Units',figUnits,...
    'Position',figPosition) ;
% loop through variables
for j = 1:length(data_mat_cell_exp)
    %--------------------
    % get current data
    data_mat_exp = data_mat_cell_exp{j} ;
    data_mat_ctrl = data_mat_cell_ctrl{j} ;
    
    %---------------------
    % take good trials
    data_mat_exp = data_mat_exp(~bad_ind_exp, :) ;
    data_mat_ctrl = data_mat_ctrl(~bad_ind_ctrl, :) ;
    
    %------------------------
    % get mean and CIs
    data_mean_exp = nanmean(data_mat_exp, 1) ;
    data_ci_exp = bootci(N_BOOT_SAMPLES,@nanmean,data_mat_exp) ;
    data_mean_ctrl = nanmean(data_mat_ctrl, 1) ;
    data_ci_ctrl = bootci(N_BOOT_SAMPLES,@nanmean,data_mat_ctrl) ;
    
    %--------------------
    % smooth results?
    if smoothFlag
        data_mean_exp = smooth(data_mean_exp) ;
        data_mean_ctrl = smooth(data_mean_ctrl) ;
        
        data_ci_exp = [smooth(data_ci_exp(1,:))' ; ...
            smooth(data_ci_exp(2,:))'] ;
        data_ci_ctrl = [smooth(data_ci_ctrl(1,:))' ; ...
            smooth(data_ci_ctrl(2,:))'] ;
    end
    %----------------
    % make subplot
    subplot(4,1,j)
    hold on
    
    %--------------------------
    % patch for opto stimulus
    ylim = ylim_cell{j} ;
    ymin = ylim(1) ;
    ymax = ylim(2) ;
    h_patch = patch(LED_start_time*[1 1 1 1] + LED_duration*[0 0 1 1],...
        [ymin ymax ymax ymin], LED_color, 'EdgeColor', 'none') ;
    h_patch.FaceAlpha = LED_alpha ;
    
    %-------------------
    % plot raw traces?
    if rawTracesFlag
        % control data
        plot(sub_time, data_mat_ctrl, '-', 'Color', [ctrlColor, line_alpha],...
            'LineWidth', lw_thin) ;
        plot(sub_time, data_mat_exp, '-', 'Color', [expColor, line_alpha],...
            'LineWidth', lw_thin) ;
    end
    
    %-------------------------
    % plot means and CIs
    
    % control
    h_fill_ctrl = fill([sub_time,fliplr(sub_time)],[data_ci_ctrl(1,:), ...
        fliplr(data_ci_ctrl(2,:))], ctrlColor,'linestyle','none');
    h_fill_ctrl.FaceAlpha = ci_alpha ;
    plot(sub_time,data_mean_ctrl,'-', 'Color',ctrlColor,...
        'LineWidth',lw_thick)
    
    % experimental
    h_fill_exp = fill([sub_time,fliplr(sub_time)],[data_ci_exp(1,:), ...
        fliplr(data_ci_exp(2,:))], expColor,'linestyle','none');
    h_fill_exp.FaceAlpha = ci_alpha ;
    plot(sub_time,data_mean_exp,'-', 'Color', expColor,...
        'LineWidth',lw_thick)
    
    %----------------------------------
    % axis properties
    set(gca,'fontSize',axisFontSize,'fontName',fontName,...
        'LineWidth',axisLineWidth)
    ylabel(ylabel_cell{j},'fontSize',labelFontSize,'fontName',fontName)
    % add xlabel if this is bottom plot
    if (j == 4)
        xlabel('Time (s)','fontSize',axisFontSize,'fontName',fontName)
    else
        set(gca,'XTickLabel',[])
    end
    
    % axis limits
    set(gca,'ylim',ylim_cell{j})
    set(gca,'xlim',xlim)
end

end

%==========================================================================
%% extra functions
function [data_mat_cell, bad_ind_cell] = fetch_data_mats(grand_mean_struct,...
    flyInd, plotVarNames_cell, excludeTrialsFlag)
%---------------------------------------------------------
% which indices to look at:
plot_ind = 75:150 ;
diff_cut_off = 2.25 ;
%---------------------------------------------------------
% get info from data structure and initialize data cell
structDataNames = {grand_mean_struct.var_name} ; % names of variables
fly_ind_mat = grand_mean_struct(1).data_ind ; % index of flies
data_mat_cell = cell(length(plotVarNames_cell),1) ;
bad_ind_cell = cell(length(plotVarNames_cell),1) ;
%------------------------------
% loop through variable types
for i = 1: length(plotVarNames_cell)
    varNamesCurr = plotVarNames_cell{i} ;
    if length(varNamesCurr) > 1
        % find structure indices for each kinematic variable type
        data_R_ind = cellfun(@(y) strcmp(y,varNamesCurr{1}),structDataNames) ;
        data_L_ind = cellfun(@(y) strcmp(y,varNamesCurr{2}),structDataNames) ;
        
        % check that we can find data
        if (sum(data_R_ind) ~= 1) || (sum(data_L_ind) ~= 1) || ...
                (sum(data_R_ind & data_L_ind) > 0)
            disp('error fetching data')
            keyboard
        end
        
        % get left and right wing data, then average
        data_R = grand_mean_struct(data_R_ind).data_mat((fly_ind_mat == flyInd),:) ;
        data_L = grand_mean_struct(data_L_ind).data_mat((fly_ind_mat == flyInd),:) ;
        
        data_mat = (data_R + data_L)./2 ;
    elseif length(varNamesCurr) == 1
        % find index for singular data type
        data_ind = cellfun(@(y) strcmp(y,varNamesCurr{1}),structDataNames) ;
        
        % check that we can find data
        if (sum(data_ind) ~= 1)
            disp('error fetching data')
            keyboard
        end
        
        data_mat = ...
            grand_mean_struct(data_ind).data_mat((fly_ind_mat == flyInd),:) ;
    else
        disp('invalid number of data types requested')
        keyboard
    end
    
    %-----------------------------------
    % exclude files that are too noisy
    if excludeTrialsFlag
        %disp('under construction')
        mean_curr = nanmean(data_mat(:,plot_ind), 1) ;
        diff_mat = bsxfun(@minus, data_mat(:,plot_ind), mean_curr) ;
        diff_mat_norm = nanmax(abs(diff_mat),[],2) ; % L1 norm
        diff_mat_zScore = (diff_mat_norm - nanmean(diff_mat_norm))./ ...
            nanstd(diff_mat_norm) ;
        bad_ind = (diff_mat_zScore > diff_cut_off) ;
        bad_ind_cell{i} = bad_ind ;
        if (0)
            figure ;
            hold on
            x = 1:length(diff_mat_zScore) ;
            plot(x,diff_mat_zScore,'ko-') ;
            plot(x(bad_ind), diff_mat_zScore(bad_ind),'rx')
            
            figure ;
            hold on
            t_curr = 1:size(data_mat,2) ;
            plot(t_curr(plot_ind), data_mat(~bad_ind, plot_ind),'b-')
            if sum(bad_ind) > 0
                plot(t_curr(plot_ind), data_mat(bad_ind, plot_ind),'r-')
            end
            axis tight
            
            %keyboard
        end
    else
        bad_ind_cell{i} = false(size(data_mat,1),1) ;
    end
    
    %-----------------------
    % store data in cell
    data_mat_cell{i} = data_mat ;
end

end