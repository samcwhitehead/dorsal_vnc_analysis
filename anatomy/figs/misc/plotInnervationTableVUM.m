% -------------------------------------------------------------------------
% script to make plot of VUM muscle innervation table
% -------------------------------------------------------------------------
%% path and params
% where to find data and save plots 
[mfilePath, ~, ~] = fileparts(mfilename('fullpath')) ; 
figDirectory = fileparts(mfilePath) ;
parentDirectory = fileparts(figDirectory) ; 
dataPath = fullfile(parentDirectory, 'data', 'cell_annotations') ;
dataFn = 'T2_VUM_Projections.xlsx' ; 

saveFlag = false ;
savePath = fullfile(mfilePath,'output') ;
saveFn = 'VUM_innervation_table.png' ;


% ----------------------------------------
% key for string to number conversion
% NB: going to use 1=innervation, 2=line has this neuron
hit_str = "+" ;   %NB: using contains, so doesn't have to be whole thing
VUM_str = "VUM" ;

% ---------------------
% plot params
gridFlag = true ; 

none_color = [1,1,1] ; 
contains_color = [0,0,0] ;
output_color = [178,24,43]/255 ; % red (color for innervation)
color_mat = [none_color ; output_color ; contains_color] ; 
    
lw_grid = 0.5 ; % line width for grid

figWidth = 3.346 ; % inches
figHeight = 7.303 ; % inches
figPosition = [7.0625, 3.7396, 1.9293, 4.3542] ; % [5.6458, 3.7396, figWidth, figHeight] ;
figUnits = 'inches' ;

fontName = 'arial' ; 
fontSize = 6 ;
fontSizeSmall = 6 ; % 4.5 ;

flipYTickLabelFlag = false ; 
XTickLabelTopFlag = true ; 
% -----------------------------------------------------------
%% read data table
% NB: a positive hit is a "+" while no innervation/expression is a ""
dataPathFull = fullfile(dataPath, dataFn) ;
try
    % detectImportOptions is causing strange errors at the moment...
    opts = detectImportOptions(dataPathFull);
    innervationTable = readtable(dataPathFull,opts) ;
catch
    innervationTable = readtable(dataPathFull) ;
end

% read out info from the table
innervCell = innervationTable{:,:} ;
targetNames = innervationTable.Properties.VariableNames(2:end) ;
lineNames = innervCell(:,1) ;

% convert cell array of "+" and "" into numbers
n_lines = length(lineNames) ; 
n_targs = length(targetNames) ; 

innervMat = zeros(n_lines, n_targs) ;
for i = 1:n_lines
    for j = 1:n_targs
        % get current target VUM or muscle
        targ_curr = targetNames{j} ;
        
        str_curr = innervCell{i, j+1} ;
        if contains(str_curr,hit_str) && contains(targ_curr, VUM_str)
            innervMat(i,j) = 2 ;
        elseif contains(str_curr,hit_str) && ~contains(targ_curr, VUM_str)
            innervMat(i,j) = 1 ;
        end
    end
end

% % flip matrix, since we have more room in height than width
%innervMat = innervMat' ;

% -------------------------------------------------------------------------
%% make plot
% initialize figure window
h_main = figure('PaperPositionMode','auto','MenuBar','none',...
    'ToolBar','none','DockControls','off','Units',figUnits,...
    'OuterPosition',figPosition) ;
ax = gca ; 
hold(ax,'on') 

% plot innervation colored squares

imagesc(ax, flipud(innervMat')) % flip matrix, since we have more room in height than width
colormap(ax, color_mat) 

% -------------------------
% plot grid lines
if gridFlag
    for i = 0:(size(innervMat,1))
        plot(ax, [i + 0.5, i + 0.5], [0.5, size(innervMat,2)+0.5], ...
            'k-', 'LineWidth',lw_grid)
    end
    for j = 0:size(innervMat,2)
        plot(ax, [0.5,size(innervMat,1)+0.5], [j + 0.5, j + 0.5 ], ...
            'k-', 'LineWidth',lw_grid)
    end
end
% ---------------------------
% set axis properties
%axis equal
xlim = [0.5, size(innervMat,1)+0.5] ; 
ylim = [0.5, size(innervMat,2)+0.5] ; 
set(ax,'xlim',xlim,'ylim',ylim)

% axis labels
set(ax,'fontName',fontName)
ax.XAxis.FontSize = fontSizeSmall ;
ax.YAxis.FontSize = fontSize ;

set(ax,'XTick',1:size(innervMat,1), ...
    'XTickLabel',lineNames,...
    'XTickLabelRotation',90)
set(ax,'YTick',1:size(innervMat,2), ...
    'YTickLabel',fliplr(targetNames))
if flipYTickLabelFlag
    set(ax,'YTickLabelRotation',180) ;
end
if XTickLabelTopFlag
    set(ax,'xaxisLocation','top')
end
pause(0.1)
try
    ax.XRuler.Axle.LineStyle = 'none' ; 
catch
    fprintf('No x axis ruler... \n')
end
try
    ax.YRuler.Axle.LineStyle = 'none' ; 
catch
    fprintf('No y axis ruler... \n')
end
    
set(ax,'TickLength',[0 0])

% ----------------------------------------------
%% save figure?
if saveFlag
   export_fig(h_main, fullfile(savePath, saveFn), '-dpng','-r600') ;  
end