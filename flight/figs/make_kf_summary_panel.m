% -------------------------------------------------------------------------
% script to create 4 axis panels for summary of kinematic DOFs from kinefly
% data
% -------------------------------------------------------------------------
%% path info and plot params
[mfilePath, ~, ~] = fileparts(mfilename('fullpath')) ; 
parentDirectory = fileparts(mfilePath) ;
dataPath = fullfile(parentDirectory, 'data') ;
savePath = fullfile(parentDirectory,'figs','output') ;

saveFlag = false ;

% order to plot different genotypes in 
driver_order = {'SS01062', 'SS44056', 'SS41068', 'SS47120', 'SS47152', ...
    'SS41039', 'SS37246', 'SS48311', 'SS37253' } ;
ctrl_driver = 'SS01062' ; 
var_type_list = {'amp', 'wbf', 'fwdDev','backDev'} ;

% --------------------------------
% plot params
figUnits = 'inches' ; 
figPosition = [3.5, 4.0, 6.375, 3.4] ; 

n_rows = 2 ; % panel layout
n_cols = 2 ; 
n_plots = n_rows*n_cols ; 

% subplot dimensions
gap = [0.02, 0.075] ; %[0.03, 0.05] ; %[0.1 0.10] ;
marg_h = [0.2, 0.01] ; %[0.08 0.08] ;
marg_w = [0.075, 0.01] ;

% ----------------------------------------------
%% load data
kf_data_struct = importdata(fullfile(dataPath, 'kf_data_struct.mat')) ;

% ------------------------------------------------
%% make figure
% initialize figure window
h_main = figure('PaperPositionMode','auto','MenuBar','none',...
    'ToolBar','none','DockControls','off','Units',figUnits,...
    'OuterPosition',figPosition) ;


% also initialize storage for axis handles
ax_array = gobjects(n_plots,1) ; 

% loop over variable types
for k = 1:length(var_type_list)
    % current kinematic dof
    var_type = var_type_list{k} ; 
    
    % row and col of current plot
    [c, r] = ind2sub([n_rows, n_cols], k) ; 
    
    % initialize axis for current dof
    ax_array(k) = subtightplot(n_rows, n_cols, k, gap, ...
        marg_h, marg_w) ;
    set(ax_array(k),'Parent', h_main)
    
    % determine whether or not we want labels on this plot
    labelFlag = (r == 2) ; 

    % make plot on current axis
    ax_array(k) = plot_kf_summary(ax_array(k), kf_data_struct, ...
        var_type, driver_order, labelFlag) ; 
end

% -----------------------------------------------------
%% save results?
if saveFlag
    % switch figure renderer
    set(h_main,'Renderer','painters') 
    
    save_fn = 'kinefly_summary_stats' ; 
    savePathFull = fullfile(savePath, save_fn) ; 
    % fig
    savefig(h_main, savePathFull) 
    % png
    exportgraphics(h_main, [savePathFull '.png'], 'Resolution', 500)
    % svg
    print(h_main, savePathFull, '-dsvg') ; 
end