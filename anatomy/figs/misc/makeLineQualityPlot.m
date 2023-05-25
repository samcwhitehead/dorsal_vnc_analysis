% -------------------------------------------------------------------------
% script to generate a plot of line quality (ranked "A" - "C") for our
% collection. NB: "A" is the best.
%
% Should brain info be included in same plot, or should we make multiple
% plots?
% -------------------------------------------------------------------------
%% path info 
[mfilePath, ~, ~] = fileparts(mfilename('fullpath')) ; 
figDirectory = fileparts(mfilePath) ;
parentDirectory = fileparts(figDirectory) ; 
dataPath = fullfile(parentDirectory, 'data') ;
dataFn = 'driver_line_list.xlsx' ; 
dataPathFull = fullfile(dataPath, dataFn) ; 

saveFlag = false ;
savePath = fullfile(mfilePath,'output') ;

% some params
ratingsVals = {'A','A/B', 'B', 'C'} ; 
bin_space = 0.5 ; % define edges as bin_center +/- "bin_space"

% plot params
bar_width = 0.6 ;
figUnits = 'inches' ; 
figPosition = [2.0313, 4.2396, 3.3229, 4.7708] ;
% ------------------------------------------------------
%% load line quality spreadsheet
try
    % detectImportOptions is causing strange errors at the moment...
    opts = detectImportOptions(dataPathFull);
    qualityTable = readtable(dataPathFull,opts) ;
catch
    qualityTable = readtable(dataPathFull) ; 
end

% doesn't seem to want to read column headers atm, so just get the
% variables that we know
lineName = qualityTable.DriverLine ; 
vncQuality = qualityTable.VNCRating ; 
brainQuality = qualityTable.BrainRating ; 

% remove repeats of headers
keep_rows = cellfun(@(y) contains(y, 'SS'), lineName) ; 
lineName = lineName(keep_rows) ; 
vncQuality = vncQuality(keep_rows) ; 
brainQuality = brainQuality(keep_rows) ; 

% -----------------------------------------------------------
%% quantify fraction of number of lines with each quality
% get indices for quality entries of each type
[vncQualUnique, ~, vncQualInd] = unique(vncQuality) ;  
[brainQualUnique, ~, brainQualInd] = unique(brainQuality) ;

% find all entries (in list of unique values) that are either "A", "B", or
% "C"
good_vnc_ind = find(cellfun(@(y) ismember(y, ratingsVals), ...
    vncQualUnique)) ; 
good_brain_ind = find(cellfun(@(y) ismember(y, ratingsVals), ...
    brainQualUnique)) ; 

% this should give us a mapping between the index (e.g. in "vncQualInd")
% and letter grade, since unique function sorts outputs. so let's define
% some edges and make this a histogram
vnc_edge_left = min(good_vnc_ind) - bin_space ; 
vnc_edge_right = max(good_vnc_ind) + bin_space ;
vnc_edges = vnc_edge_left:vnc_edge_right ; 

brain_edge_left = min(good_brain_ind) - bin_space ; 
brain_edge_right = max(good_brain_ind) + bin_space  ; 
brain_edges = brain_edge_left:brain_edge_right ; 

% get hist counts
vnc_qual_counts = histcounts(vncQualInd, vnc_edges) ; 
brain_qual_counts = histcounts(brainQualInd, brain_edges) ; 

% ---------------------------------------------------------------
%% go through and get joint quality ratings
dataMat = zeros(length(good_vnc_ind), length(good_brain_ind)) ; 
for k = 1:length(good_vnc_ind)
    % index for lines with current quality rating
    idx = (vncQualInd == good_vnc_ind(k)) ; 
    
    % get their count of brain ratings
    dataMat(k,:) = histcounts(brainQualInd(idx), brain_edges) ; 
end
% ---------------------------------------------------------------
%% make plots
% NB: it's probably more informative to break up the vnc bars by brain
% quality, as opposed to doing three different plots
h_main = figure('PaperPositionMode', 'auto', 'Units', figUnits,...
    'OuterPosition',figPosition) ;
ax = gca ; 
hold on

% make bar plot
b = bar(ax, dataMat, 'stacked') ; 

% adjust bar properties (can add more later, but for now just removing
% baseline)
for b_num = 1:length(b)
    % remove baseline
    b(b_num).BaseLine.LineStyle = 'none';
    
    % set bar width
    b(b_num).BarWidth = bar_width ; 
end

% --------------------------------------
% axis properties
% --------------------------------------
% set axis ticks and labels
x_ticks = 1:length(good_vnc_ind) ; 
set(ax, 'XTick', x_ticks, 'XTickLabel', ratingsVals) 
set(ax, 'xlim', [x_ticks(1) - bar_width/2, x_ticks(end) + bar_width/2]) 
xlabel('VNC sparsity rating')
ylabel('driver line count')


% add in legend
lgd = legend({'A', 'B', 'C'}, 'location', 'northwest') ;
title(lgd,'brain sparsity rating', 'FontWeight','normal') 
legend('boxoff')

% % make axis look similar to other plots
% ax = b1_paper_prettify_axis(ax, true, false, false, false, true, ...
%     [1.25, 0.85]) ; 

% -----------------------------------------------------------------
%% save output?
if saveFlag
    savePathFull = fullfile(savePath, 'line_quality_hist') ; 
    
    savefig(h_main, [savePathFull '.fig']) ; 
    exportgraphics(h_main, [savePathFull '.png'], 'Resolution', 500) ; 
    print(h_main, [savePathFull '.svg'], '-dsvg') ; 
end



