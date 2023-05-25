% -------------------------------------------------------------------------
% function to plot kf summary data to an existing axis
% -------------------------------------------------------------------------
function ax = plot_kf_summary(ax, kf_data_struct, var_type, ...
    driver_order, labelFlag)
% ---------------------
%% inputs and params
if ~exist('ax','var') || isempty(ax)
    ax = gca ;
end
if ~exist('kf_data_struct','var') || isempty(kf_data_struct)
    [mfilePath, ~, ~] = fileparts(mfilename('fullpath')) ; 
    parentDirectory = fileparts(mfilePath) ;
    dataPath = fullfile(parentDirectory, 'data') ; 
    kf_data_struct = importdata(fullfile(dataPath, 'kf_data_struct.mat')) ;
end
if ~exist('var_type','var') || isempty(var_type)
    var_type = 'amp' ; % 'amp' | 'fwdDev' | 'backDev' | 'wbf' ;
end
if ~exist('driver_order','var') || isempty(driver_order)
    driver_order = {'SS01062', 'SS44056', 'SS41068', 'SS47120', 'SS47152', ...
        'SS41039', 'SS37246', 'SS48311', 'SS37253' } ;
end
if ~exist('labelFlag','var') || isempty(labelFlag)
    labelFlag = true ;
end

% name of empty split being used for control group
ctrl_driver = 'SS01062' ;

% remove outliers?
removeOutliersFlag = false ; 

% % show multcompare stats?
% displayStatsStr = 'on' ; % 'off' | 'on' ; 

% delta character 
delta_char = char(hex2dec('0394')) ;
% ---------------------------------
% general plot params
define_constants_kf

% % bar plots -- i guess just using defaults?
% bar_width = 0.7 ; % 0.5
% barEdgeColor = 'k' ; 
% barAlphaVal = 0.33 ;

% zero line
lineWidthThin = 0.75 ; 
lineStyle = '--' ; 
lineColor = 'k' ; 

% data plots
showMM = 'median' ; 

% x axis ticks and limits (y axis depends on data type)
x_ticks = 1:length(driver_order) ;
xlim = [0.5, length(x_ticks) + 0.5] ;

% axis hold situation
holdFlag = ishold(ax) ;
hold(ax,'on') ;

% ----------------------------------------------------------
%% select kinematic DOF
% get field names (as well as other associated info) for this DOF
switch var_type
    case 'amp'
        field_names = {'R_AMP', 'L_AMP'} ;
        ylim = 12*[-1, 1] ;
        dof_unit = 'deg' ;
        dof_name = 'stroke amp' ; 
        
    case 'fwdDev'
        field_names = {'R_DEV_1', 'L_DEV_2'} ;
        ylim = 7*[-1, 1] ;
        dof_unit = 'deg' ;
        dof_name = 'fwd deviation' ; 
        
    case 'backDev'
        field_names = {'R_DEV_2', 'L_DEV_1'} ;
        ylim = 7*[-1, 1] ;
        dof_unit = 'deg' ;
        dof_name = 'back deviation' ; 
        
    case 'wbf'
        field_names = {'WBF'} ;
        ylim = 12*[-1, 1] ;
        dof_unit = 'Hz' ;
        dof_name = 'wingbeat freq' ; 
end

% ----------------------------------------------------------------------
%% grab data from struct
% initialize storage for data and group info
data = [] ;
grp_ind = [] ;
grp_names = [] ;

% clumsily loop through to grab all data (need in vector form for stats)
for k = 1:length(driver_order)
    % index in struct for current driver
    driver_curr = driver_order{k} ;
    driver_idx = arrayfun(@(x) strcmpi(x.driver, driver_curr), ...
        kf_data_struct) ;
    
    % check if driver is there
    if sum(driver_idx) ~= 1
        fprintf('Could not find data for %s \n', driver_curr)
        keyboard
    end
    
    % --------------------------------------
    % read out data for current driver
    for m = 1:length(field_names)
        % current field name
        f_name = [field_names{m} '_peak_mat'] ;
        
        % compile current data
        if m == 1
            data_curr = kf_data_struct(driver_idx).(f_name) ;
        else
            data_curr = data_curr + kf_data_struct(driver_idx).(f_name) ;
        end
    end
    
    % divide by 2 to make average (if data is L/R)
    data_curr = data_curr./length(field_names) ;
    
    % add data to compiled vector
    data = [data ; data_curr] ;
    
    % ---------------------------------------
    % get group info
    grp_ind_curr = k.*ones(size(data_curr)) ;  % integer index
    grp_names_curr = cell(size(data_curr)) ;  % driver name grouping
    grp_names_curr(:) = {driver_curr} ;
    
    % add group info to compiled arrays
    grp_ind = [grp_ind ; grp_ind_curr] ;
    grp_names = [grp_names ; grp_names_curr] ;
    
end

% ---------------------------------------------------------------
%% get control data for stats
ctrl_driver_ind = find(cellfun(@(y) strcmpi(ctrl_driver, y), driver_order));
ctrl_data = data(grp_ind == ctrl_driver_ind) ; 

% -----------------------------------------------------------------------
%{
% *** Switched from multcompare because it made k*(k-1)/2 comparisons
% when I only need (k-1) for each experimental group vs control. Now using
% pairwise ranksum with bonferroni correction

%% run stats for current data
% % run kruskal wallis test
% [p, ~, stats] = kruskalwallis(data, grp_names, displayStatsStr) ;
% fprintf('Kruskal wallis for %s: p = %f \n', var_type, p)
% 
% % run multcompare
% mc_mat = multcompare(stats,'CType','tukey-kramer','Display',displayStatsStr) ;
% 
% % get index in data struct for controller struct
% ctrl_driver_ind = find(cellfun(@(y) strcmpi(ctrl_driver, y), driver_order));

% % get pval list (one entry per driver)
% pValList = ones(size(driver_order)) ; 
% for q = 2:length(driver_order)
%    % find entry for current driver in matrix "mc_mat" (from multcompare)
%    row_idx = (mc_mat(:,1) == ctrl_driver_ind) & (mc_mat(:,2) == q) ; 
%    if sum(row_idx) ~= 1
%       fprintf('Error: could not find p value for driver: %s \n', ...
%           driver_order{q}) 
%       keyboard
%    end
%    
%    % add current p value to list
%    pValList(q) = mc_mat(row_idx,end) ; 
% end
%}
% -------------------------------------------------------------------
%% plot data to axis 
% plot line at y = 0
plot(ax, xlim, [0,0], lineStyle,...
    'Color', lineColor, ...
    'LineWidth', lineWidthThin) ; 

% loop over drivers
for j = 1:length(driver_order)
    % current driver
   driver_curr = driver_order{j} ;  
   
   % get associated data
   data_idx = (grp_ind == j) ; 
   data_curr = data(data_idx) ; 
   
   % remove outliers?
   if removeOutliersFlag
       data_curr = rmoutliers(data_curr) ; 
   end
   
   % also get associated p value
   % ***NB: switched from multcompare (which made k*(k-1)/2 comparisons
   % when I only need (k-1) for each experimental group vs control. so now
   % using pairwise ranksum with bonferroni correction
   pVal = ranksum(data_curr, ctrl_data) ; % get wilcoxon p value
   pVal = pVal*(length(driver_order) - 1) ; % apply bonferroni correction
   
   % change coloring and filling based on whether or not it's
   % significant
   if (pVal <= 0.05) && ~strcmp(driver_curr, ctrl_driver)
       cmap_curr = brewermap(6,cmap_struct.SS37246) ;
       plotColor = cmap_curr(5,:) ;
   elseif (pVal > 0.05) && ~strcmp(driver_curr, ctrl_driver)
       plotColor = 'k' ;
   else
       cmap_curr = brewermap(6,'Greys') ;
       plotColor = cmap_curr(5,:) ;
   end
                    
   % --------------------------------------------
   % plot data
   ax = myPlotGroupDataToAxis(ax, j, data_curr, ...
       'showMM', showMM, ...
       'plotColor', plotColor, ...
       'barColor', plotColor) ;
   
   % add significance star(s)?
   if pVal < 0.001
       text(j, ylim(2), '***','Units','data', ...
           'HorizontalAlignment', 'center')
   elseif pVal < 0.01
       text(j, ylim(2), '**','Units','data', ...
           'HorizontalAlignment', 'center')
   elseif pVal < 0.05
       text(j, ylim(2), '*','Units','data', ...
           'HorizontalAlignment', 'center')
   end
   
   % add label for driver?
   if labelFlag
       mn_curr = mn_name_struct.(driver_curr) ;
%        text_str = sprintf('%s\n({{\\bf%s}})',driver_curr, mn_curr) ;
       text_str = sprintf('%s ({{\\bf%s}})',driver_curr, mn_curr) ;
       text(j , ylim(1) - 0.65, text_str, 'fontName',fontName,...
           'fontSize', labelFontSize, ...
           'fontWeight','normal', ...
           'HorizontalAlignment', 'right',...
           'VerticalAlignment', 'middle',...
           'Rotation', 45) ;
   end
end

% ------------------------------------------------------------
%% adjust axis properties
% set axis font and limit properties
set(ax,'fontName',fontName,'fontSize',axisFontSize)
set(ax,'xlim',xlim,'ylim',ylim)
set(ax,'Color','none') 

% set x ticks (in case we decide to label them)
set(ax, 'XTick', x_ticks)

% add y label 
ylabel(ax, sprintf('%c %s\n(%s)', delta_char, dof_name, dof_unit))

% make axis look nice
ax = b1_paper_prettify_axis(ax, true, false, true, false, true, 0.75) ; 

% adjust x axis ruler (depends on how we want to label)
set(ax,'xcolor','none')

% return axis to previous hold state
if ~holdFlag
   hold(ax, 'off') 
end

end