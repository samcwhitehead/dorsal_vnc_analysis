% -------------------------------------------------------------------------
% function to plot a vnc surface to the current axis -- to be used with
% neuron plots
% -------------------------------------------------------------------------
function ax = drawVNCBackgroundToAxis(ax, fig, vnc_struct, vnc_color, ...
    vnc_alpha, parent)
% --------------------
%% inputs and params
if ~exist('ax','var') || isempty(ax)
    ax = gca ; 
end
if ~exist('fig','var') || isempty(fig)
    fig = gcf ; 
end
if ~exist('vnc_struct','var') || isempty(vnc_struct)
    dataPath = 'D:\Fly Imaging\Erica Interneuron Stacks\VNC_for_drawings\' ;
    dataName = 'vnc_struct.mat' ; 
    vnc_struct = importdata(fullfile(dataPath, dataName)) ; 
end
if ~exist('vnc_color','var') || isempty(vnc_color) 
    vnc_color = 0.7*[1, 1, 1] ;
end
if ~exist('vnc_alpha','var') || isempty(vnc_alpha) 
    vnc_alpha = 0.2 ;
end
if ~exist('parent', 'var') || isempty(parent)
    parent = [] ; 
end

% check hold status of current axis
set(fig,'CurrentAxes',ax)
holdFlag = ishold(ax) ; 
hold(ax,'on') ; 

% -----------------------------------
%% read data from vnc struct and plot
% points, boundary, and shrink factor
qx = vnc_struct.qx ; 
qy = vnc_struct.qy ;
qz = vnc_struct.qz ; 
vnc_bound = vnc_struct.vnc_bound ; 
sf = vnc_struct.shrinkFactor ; 

% make triangulation surface plot -- scale back up by scale factor
vncSurf = trisurf(vnc_bound, qx*sf, qy*sf, qz*sf, ...
    'FaceColor',vnc_color, 'EdgeColor','none','FaceAlpha',vnc_alpha) ; 

% --------------------------------------
%% axis properties
axis(ax, 'equal')
if holdFlag
    hold(ax ,'on')
else
    hold(ax, 'off')
end

% parent properties
if ~isempty(parent)
    grp = hgtransform('Parent',parent) ;
    set(vncSurf,'Parent',grp) ;  
end

end

