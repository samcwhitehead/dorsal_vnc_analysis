% -------------------------------------------------------------------------
% function to generate 3D plots of segmented interneurons as sanity check 
% on overlap analysis. Assumes images are binarized and the same size
% -------------------------------------------------------------------------
function h_IN = plotMultipleNeurons(bwCell, labels)
% -----------------
% inputs
if ~exist('labels','var') 
   labels = [] ;  
end
% -----------------
% plot params
N_images = length(bwCell) ; 
color_mat = lines(N_images) ; 

marker_size = 6 ; 
plot_alpha = 0.2 ; 

% ---------------
% make plot
h_IN = figure('PaperPositionMode','auto') ;
hold on
ax_array = gobjects(N_images) ;

for i = 1:N_images
    bwMat = bwCell{i} ;
    [x, y, z] = ind2sub(size(bwMat),find(bwMat(:))) ;
    %plot3(x, y, z, '.','Color',[color_mat(i,:), 0.01])
    ax_array(i) = plot3(nan, nan, nan, '.','MarkerSize', 3*marker_size,...
        'Color',color_mat(i,:)) ;
    h = scatter3(x, y, z, marker_size, color_mat(i,:), '.','HandleVisibility','off') ;
    set(h, 'MarkerEdgeAlpha', plot_alpha, 'MarkerFaceAlpha', plot_alpha)
    
end

axis equal
axis off

if ~isempty(labels)
   [l, icons] = legend(labels, 'location','northwest','Interpreter','none') ;
   
end

end