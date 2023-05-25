% -------------------------------------------------------------------------
% should take in the output of "mySkelAndGraph.m" and plot results in a way
% that looks reasonable
% -------------------------------------------------------------------------
function ax = drawSkelGraph(ax, skel_struct, lineColor, lineWidth, ...
    lineAlpha, smoothFlag, flipYFlag)
% inputs and plot params
if ~exist('ax','var') || isempty(ax)
    ax = gca ;
end
if ~exist('lineColor','var') || isempty(lineColor)
    lineColor = 0.5*[1,1,1] ;
end
if ~exist('lineWidth','var') || isempty(lineWidth)
    lineWidth = 0.5 ;
end
if ~exist('lineAlpha','var') || isempty(lineAlpha)
    lineAlpha = 1.0 ;
end
if ~exist('smoothFlag','var') || isempty(smoothFlag)
    smoothFlag = false ;
end
% normally we want this (image coordinate correction), but when drawing
% alongside VNC this messes things up if true
if ~exist('flipYFlag','var') || isempty(flipYFlag)
    flipYFlag = true ;
end

% -----------------------------------------
% make sure axis is set to "hold on"
holdFlag = ishold(ax) ; 
hold(ax,'on') ; 

% -----------------------------------------
% read data
% image dimensions
[W,L,H] = size(skel_struct.skel) ;

% read out links and nodes
links = skel_struct.link ;
nodes = skel_struct.node ;

% -------------------------------------------------------------------------
% loop over links and draw each. include the node endpoints, and use an
%  interpolant to smooth the lines through them
for i = 1:length(links)
    % read out link points
    [x, y, z] = ind2sub([W,L,H], links(i).point) ;
    
    % append node endpoints
    n1 = links(i).n1 ;
    n2 = links(i).n2 ;
    x_all = [nodes(n1).comx, x, nodes(n2).comx] ;
    y_all = [nodes(n1).comy, y, nodes(n2).comy] ;
    z_all = [nodes(n1).comz, z, nodes(n2).comz] ;
    
    % flip for image coordinates?
    if flipYFlag
        x_all = L - x_all ; % this is labeled x but will become y when plotting :(
    end
    
    % get points in single array
    pts_all = [y_all', x_all', z_all'] ;
    
    %    % spline through these points
    %    curve = cscvn(pts_all') ;
    %    curve_pts = fnval(pts_all, curve) ;
    
    % smooth links?
    if smoothFlag
        for dim = 1:size(pts_all,2)
            pts_all(:,3) = smooth(pts_all(:,3)) ;
        end
    end
    % plot link
    plot3(pts_all(:,1), pts_all(:,2), pts_all(:,3), '-',...
        'Color', [lineColor, lineAlpha], 'LineWidth', lineWidth)
    
    
end

% set some axis properties
axis image;axis off;
set(gcf,'Color','w');
drawnow;

% return axis to previous hold state
if holdFlag
    hold(ax ,'on')
else
    hold(ax, 'off')
end

end