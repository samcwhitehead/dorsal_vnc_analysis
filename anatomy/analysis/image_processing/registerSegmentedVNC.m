% -------------------------------------------------------------------------
% script to register the segmented VNC + subregion images that Hiro sent me
% to the standard 2018 VNC Janelia template coordinates I've been using
% -------------------------------------------------------------------------
%% get path info
dataRoot = 'D:\Fly Imaging\Erica Interneuron Stacks\' ;
segVNCPath = fullfile(dataRoot , 'VNC neuropil', 'segmented_images') ; % oath to images that i segmented
tformPath = fullfile(segVNCPath, 'registered') ;
tformFn = 'tform_cell.mat' ;

savePath = tformPath ;

overWriteFlag = true ;
% plotFlag = true ;

boundaryConvexity = 1.0 ; % 0 gives convex hull, 1 compact boundary

% get directory of segmented VNC regions
segVNCDir = dir(fullfile(segVNCPath, '*.mat')) ;

% load tform cell
tform_cell = importdata(fullfile(tformPath, tformFn)) ;

% ----------------------------------------------------
%% loop over neuropil regions and apply transforms
for k = 1:length(segVNCDir)
    % ---------------------------------------
    % load current neuropil region data
    pathCurr = fullfile(segVNCDir(k).folder, segVNCDir(k).name) ;
    neuropil_curr = importdata(pathCurr) ;
    [~, neuropil_name, ~] = fileparts(pathCurr) ;
    
    % number of components in neuropil image
    N_objs = length(neuropil_curr.vox_coord_cell) ;
    % initialize output storage
    vox_coord_cell = cell(N_objs,1) ;
    vox_boundary_cell = cell(N_objs,1) ;
    
    % loop over components of neuropil image
    for m = 1:N_objs
        % convert coordinates to point cloud to make tform easier
        pc_np = pointCloud(neuropil_curr.vox_coord_cell{m}) ;
        
        % ----------------
        % apply tforms
        for n = 1:length(tform_cell)
            % get current tform
            tform = tform_cell{n} ;
            
            % apply transformation
            pc_np = pctransform(pc_np, tform) ;
        end
        
        % ---------------------------------------------------------------------
        % recalculate boundary and add current data to output cells
        % coordinates
        vox_coord_cell{m} = pc_np.Location ;
        
        % re-calculate boundary
        vox_boundary = boundary(pc_np.Location(:,1), pc_np.Location(:,2),...
            pc_np.Location(:,3), boundaryConvexity) ;
        
        vox_boundary_cell{m} = vox_boundary ;
    end
    
    % ---------------------------------------------------------
    % save registered neuropil data
    savePathCurr = fullfile(savePath, [neuropil_name '.mat']) ;
    if ~exist(savePathCurr,'file') || overWriteFlag
        save(savePathCurr, 'vox_coord_cell', 'vox_boundary_cell')
    end
    
    fprintf('Completed %d/%d neuropil regions \n', k, length(segVNCDir))
end
