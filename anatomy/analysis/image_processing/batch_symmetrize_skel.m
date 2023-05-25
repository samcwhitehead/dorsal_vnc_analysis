% -------------------------------------------------------------------------
% script to go through non-symmetrized skeleton structures and generate
% their symmetrized counterparts
% -------------------------------------------------------------------------
%% data path info
dataRoot = 'D:\Fly Imaging\Erica Interneuron Stacks\IN\' ;
dataPath = fullfile(dataRoot, 'skeletonized_non_sym') ; 
savePath = fullfile(dataRoot, 'skeletonized_sym') ; 
if ~exist(savePath,'dir')
    mkdir(savePath)
end

skelDir = dir(fullfile(dataPath, '*bw_skel_struct.mat')) ;
N_skel = length(skelDir) ;

% save over? plot?
overWriteFlag = false ;
plotFlag = true ;
debugFlag = false ; 

% threshold for determining whether we have a unilateral image or not
symThresh = 0.3 ; % 0.3 for interneurons, 0.05 for DNs

% ----------------------------------------------
% loop over images and skeletonize them
for ind = 1:N_skel
    tic
    
    fprintf('Processing %d / %d stacks \n', ind, N_skel)
    % get filename for image
    dataFilename = fullfile(skelDir(ind).folder, skelDir(ind).name) ;
    [~, fn, ext] = fileparts(dataFilename) ;
    % set save path
    skeletonSavePath = fullfile(savePath, [fn '_skel_struct.mat']) ;
    
    % check if we've already converted this one
    if (~overWriteFlag) && exist(skeletonSavePath,'file')
        fprintf('Already analyzed file: %s \n', dataFilename)
        continue
    end
    % -----------------------------------------------------------
    % load skeleton struct
    skel_struct = importdata(dataFilename) ;
    
    % ---------------------------------------------------------------------
    % check if image is already bilateral
    skelMat = skel_struct.skel ; 
    [W,L,H] = size(skelMat) ; % size of skeleton image
    symCheck = sum(skelMat & flip(skelMat,2),'all')/sum(skelMat(:)) ;
    
    % if it passes the symmetry check, no need to symmetrize
    if symCheck > symThresh
        fprintf('Image %s already symmetrized \n', fn)
        continue
    end
    
    % ---------------------------------------------------------------------
    % if we get here, skeleton needs to be symmetrized. first find midpoint
    % of image on "x axis," then loop through nodes and links
    midline_x = round(L/2) ; 
    
    % initialize new data arrays for nodes, links, and skeleton
    node_new = skel_struct.node ; 
    link_new = skel_struct.link ;
    % NB: flipping skeleton here, so no need to process elsewhere
    skel_new = flip(skel_struct.skel,2) ; 
    
    % loop over nodes
    for i=1:length(node_new)
        % first reflect node points
        node_new(i).comx = node_new(i).comx;
        node_new(i).comy = -1*(node_new(i).comy - midline_x) + midline_x ;
        node_new(i).comz = node_new(i).comz;
    end    
    % loop over links
    for j=1:length(link_new)
        
        % get points associated with current link
        points_curr = [link_new(j).point] ;
        
        % convert from index to sub, do refection, then convert back to
        % index
        [x,y,z] = ind2sub([W,L,H], points_curr) ;
        x_new = x ;
        y_new = -1*(y - midline_x) + midline_x ;
        z_new = z ;
        points_new = sub2ind([W,L,H], x_new, y_new, z_new) ;
        
        
        % assign new points to link structure
        link_new(j).point = points_new ;
        
%         % also switch endpoint node indices
%         temp = link_new(j).n1 ;
%         link_new(j).n1 = link_new(j).n2 ;
%         link_new(j).n2 = temp ;
        
        
    end

    % -------------------------------------------------------------------
    % add new skeleton to original struct
    skel_struct_old = skel_struct ; 
    new_ind = length(skel_struct_old) + 1 ; 
    
    skel_struct(new_ind).node = node_new ; 
    skel_struct(new_ind).link = link_new ; 
    skel_struct(new_ind).skel = skelMat ; 
    
    % --------------------------------------------------
    % debug results?
    if debugFlag || plotFlag
        % plot old and new skeleton
        h_temp = figure ; 
        ax = gca ;
        hold on
        ax = drawSkelGraph(ax,skel_struct(1), [0,0,1]) ;
        ax = drawSkelGraph(ax,skel_struct(2), [1,0,0]) ;
        
        % make sum dummy plot objects to have legend
        h_old = plot3(nan, nan, nan, 'b-') ;
        h_new = plot3(nan, nan, nan, 'r-') ;
        legend([h_old, h_new], {'old', 'new'}, 'location','northwest')
        
        % add title
        title(fn,'Interpreter','none') 
        
        % either save or pause, depending on flags
        if debugFlag
            keyboard
        end
        
        if plotFlag 
            saveFigPath = fullfile(savePath, [fn '_skel_plot.png']) ;
            print(h_temp, saveFigPath, '-dpng', '-r300')
        end
        
        % close figure
        close(h_temp) 
        
    end
    
    % -------------------------------------------------------------
    % save resulting structure
    save(skeletonSavePath, 'skel_struct') 
    toc
end