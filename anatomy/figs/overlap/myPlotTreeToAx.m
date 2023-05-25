% -------------------------------------------------------------------------
% function to make a plot that shows connections between a given DN (or set
% of DNs) and the INs they are "connected to" via volume overlap
%
% INPUTS:
%   - ax: axis object to plot to
%   - parentNodes: string giving the name of the parent DN
%   - childNodes: Nx1 cell array containing names of children (INs)
%   - connections: 1xN matrix containing connection "stengths" (overlap)
%   - maxNumNodes: maximum number of children nodes (if input is greater,
%       reduce to strongest connections)
%   - edgeSortFlag: boolean. sort child nodes by edge strength?
%   - cmap_struct: struct containing info to color edges inidividually. if
%       left empty, defaults to "arrowColor"
%   - nodeColor: rgb vector containing node color
%   - arrowColor: rgb vector containing arrow color
%   - arrowAlpha: double value in [0,1] giving alpha val for edges
%   - fontSize: size for labels on nodes
%
% OUTPUTS:
%   - ax: axis object containing tree plot
%   - hg: graph plot object
% -------------------------------------------------------------------------
function [ax, hg] = myPlotTreeToAx(ax, parentNode, childNodes, ...
    connections, maxNumNodes, edgeSortFlag, cmap_struct, nodeColor, ...
    arrowColor, arrowAlpha,  fontSize)
% -----------------------------
%% inputs and params
if ~exist('ax','var') || isempty(ax)
    ax = gca ;
end
if ~exist('maxNumNodes','var') || isempty(maxNumNodes)
    maxNumNodes = 20 ;
end
if ~exist('edgeSortFlag','var') || isempty(edgeSortFlag)
    edgeSortFlag = true ;
end
if ~exist('cmap_struct','var') || isempty(cmap_struct)
    cmap_struct = [] ;
end
if ~exist('nodeColor','var') || isempty(nodeColor)
    nodeColor = [0, 0, 0] ;
end
if ~exist('arrowColor','var') || isempty(arrowColor)
    arrowColor = [0, 0, 0] ;
end
if ~exist('arrowAlpha','var') || isempty(arrowAlpha)
    arrowAlpha = 0.8 ; % 0.5
end
if ~exist('fontSize','var') || isempty(fontSize)
    fontSize = 6 ;
end

% other plot params
fontName = 'arial' ;
nodeMarkerSize = 4 ;
nodeMarkerType = 'o' ;
maxArrowWidth = 2 ;
minArrowWidth = 0.25 ;
arrowPosition = 0.85 ;
arrowSize = 7 ;

textSpace = 0.1 ;

y_val_parent = 1 ;
y_val_children = 0 ;
x_range = [0, 1] ;

addEllipsisFlag = false ;

% set axis hold properties
axes(ax) ;
wasOnHold = ishold ;
hold on ;

% ----------------------------------------------------
%% normalize connections and prune any that are empty
% prune empty connections
empty_idx = (connections < 1) ;
connections = connections(~empty_idx) ;
childNodes = childNodes(~empty_idx) ;

% convert connection strength (voxel count) to arrow width
max_connection = max(connections) ;
min_connection = min(connections) ;

edgeWidths = (connections - min_connection)./max_connection ;
edgeWidths = maxArrowWidth.*edgeWidths + minArrowWidth ;

% ---------------------------------------------------------
%% restrict number of child nodes, if necessary 
if length(childNodes) > maxNumNodes
    % take only the strongest connections
    [~, sort_ind] = sort(connections, 'descend') ; 
    keep_idx = false(size(connections)) ; 
    keep_idx(sort_ind(1:maxNumNodes)) = true ; 
    childNodes = childNodes(keep_idx) ; 
    connections = connections(keep_idx) ; 
    edgeWidths = edgeWidths(keep_idx) ; 
    
    % add ellipses to plot to indicate that we are leaving some out
    addEllipsisFlag = true ;
end

% ---------------------------------------------------------
%% sort child nodes by edge size?
if edgeSortFlag
    [~, sort_ind] = sort(connections, 'descend') ;
    childNodes = childNodes(sort_ind) ;
    connections = connections(sort_ind) ;
    edgeWidths = edgeWidths(sort_ind) ;
end
% ---------------------------------------------------------
%% get x positions of nodes
% first get x positions of child and parent nodes
x_val_parent = mean(x_range) ;
N_children = length(childNodes) ;
x_val_children = linspace(x_range(1), x_range(2), N_children) ;

% -----------------------------------------------------------
%% create directed graph from parent and child nodes
% generate adjaceny matrix from "connections" vector
% NB: Node 1 will be the parent; Nodes 2-(N+1) will be the children
A = zeros(N_children + 1) ;
A(1,2:end) = edgeWidths ;

G = digraph(A) ;

% ---------------------------------------------------
%% plot directed graph
% NB: if we want to use a colormap, need to use EdgeCData instead of
% EdgeColor
if ~isempty(cmap_struct)
    % in this case, we want to color each edge differently (according to
    % cmap)
    hg = plot(ax, G, 'LineWidth',edgeWidths, ...
        'Xdata', [x_val_parent, x_val_children],...
        'Ydata', [y_val_parent, repmat(y_val_children,1, N_children)]) ;
    
    % read out cmap info
    cmap = cmap_struct.cmap ; 
    clim = cmap_struct.clim ; 
    clim_edges = linspace(clim(1), clim(end), size(cmap,1)) ; 
    
    % loop over edges to assign colors
    for k = 1:N_children
        % get current connection 
        conn_curr = connections(k) ; 
        
        % get color based on connection (overlap)
        if conn_curr >= clim(end)
            color_curr = cmap(end,:) ; 
        else
            [~, ~, bin] = histcounts(conn_curr, clim_edges) ; 
            color_curr = cmap(bin,:) ; 
        end
        
        % highlight edge to apply color
        highlight(hg, [1, (k+1)], 'EdgeColor', color_curr)
    end
    
    % since this seems to mess up svg saving of arrows, remove them 
    hg.ShowArrows = false ; 
    
else
    % in this case, all edges get the same color (arrowColor)
    hg = plot(ax, G, 'LineWidth',edgeWidths, ...
        'EdgeColor', arrowColor, ...
        'Xdata', [x_val_parent, x_val_children],...
        'Ydata', [y_val_parent, repmat(y_val_children,1, N_children)]) ;
end

% apply general settings
hg.NodeColor = nodeColor ; 
hg.Marker = nodeMarkerType ; 
hg.MarkerSize = nodeMarkerSize ; 
hg.NodeLabel = [] ; 
hg.ArrowPosition = arrowPosition ; 
hg.ArrowSize = arrowSize ; 
hg.EdgeAlpha = arrowAlpha ; 


% ----------------------------------------------------------
%% add node labels
% first add label to parent node
parent_label = text(ax, x_val_parent, y_val_parent + textSpace, parentNode, ...
    'fontSize', fontSize, ...
    'fontName', fontName, ...
    'HorizontalAlignment', 'center',...
    'VerticalAlignment', 'bottom') ;

% then add children labels
children_labels = gobjects(N_children,1) ;
for k = 1:N_children
    % make text object
    cl = text(ax, x_val_children(k), y_val_children - textSpace,...
        childNodes{k},...
        'fontSize', fontSize, ...
        'fontName', fontName, ...
        'Rotation', 90, ...
        'HorizontalAlignment', 'right',...
        'VerticalAlignment', 'middle') ;
    
    % store text object for later
    children_labels(k) = cl ;
end

% add ellipsis if we had to cut out some neurons
if addEllipsisFlag
    ellipsis_x = x_val_children(end) + mean(diff(x_val_children))/2 ; 
   et = text(ax, ellipsis_x , y_val_children, '...',...
        'fontSize', 8, ...
        'fontName', fontName, ...
        'HorizontalAlignment', 'left',...
        'VerticalAlignment', 'baseline') ;
end
% --------------------------------------------------
%% set some axis properties
% remove axis ruler, ticks, and labels
axis off

% revert to original hold state
if (~wasOnHold)
    hold off ;
end
end