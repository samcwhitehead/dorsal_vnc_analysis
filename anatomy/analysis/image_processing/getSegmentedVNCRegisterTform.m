% -------------------------------------------------------------------------
% script to register the segmented VNC + subregion images that Hiro sent me
% to the standard 2018 VNC Janelia template coordinates I've been using
% -------------------------------------------------------------------------
%% get path info squared away
dataRoot = 'D:\Fly Imaging\Erica Interneuron Stacks\' ;
rawVNCPath = fullfile(dataRoot, 'VNC neuropil') ; % path to tif of segmented vnc
segVNCPath = fullfile(rawVNCPath , 'segmented_images') ; % oath to images that i segmented
% rawVNCFn = 'confocal_binary.tif' ;
segVNCFn = 'VNC.mat' ;

refVNCPath = fullfile(dataRoot, 'VNC_for_drawings') ; % path to example 2018 Janelia template vnc
refVNCFn = 'vnc_struct.mat' ;

savePath = fullfile(segVNCPath, 'registered') ;
if ~exist(savePath,'dir')
    mkdir(savePath)
end

grid_ss_seg = 10.0 ; % how much to subsample segmented image
grid_ss_ref = 2.0 ; % how much to subsample reference image
debugFlag = false ;
plotFlag = true ;
overWriteFlag = false ;
% -------------------------------------------
%% load data
fprintf('Loading and processing images...\n')

% load VNC from Hiro (pre-processed)
segVNC_struct = importdata(fullfile(segVNCPath, segVNCFn)) ;
coords_seg = segVNC_struct.vox_coord_cell{1} ;


% load VNC to register to
refVNC_struct = importdata(fullfile(refVNCPath, refVNCFn)) ;
coords_ref = (refVNC_struct.shrinkFactor).*[refVNC_struct.qx, ...
    refVNC_struct.qy, refVNC_struct.qz] ;


% convert both sets of coordinates to point clouds
pc_seg = pointCloud(coords_seg) ;
pc_ref = pointCloud(coords_ref) ;

clear coords_ref coords_seg
% -------------------------------------------------------------------------
%% try to make an initial guess at the  transformation (seg -> ref)
fprintf('Estimating initial transformation ... \n')

% initialize storage for transformations
tform_cell = {} ;

% ------------------------------------------------------
% get translation vector needed to move seg to origin
t_vec = -1.*mean(pc_seg.Location) ;

% turn translation vector into tform
tform_t1 = affine3d([1, 0, 0, 0; 0, 1, 0, 0 ; 0, 0, 1, 0 ; ...
    t_vec, 1]) ;

% apply translation
pc_seg_new = pctransform(pc_seg, tform_t1) ;

% store tform
tform_cell{end+1} = tform_t1 ;

% if debugFlag
%     figure ;
%     pcshowpair(pc_ref, pc_seg_new)
%     keyboard
% end
% ---------------------------------------------------------
% rotate seg so its long axis aligns with y axis
coefs = pca(pc_seg_new.Location) ;
phi = atan2(coefs(2,1), coefs(1,1)) - 3*pi/2 ;
R = eulerRotationMatrix(-phi, 0, 0) ;

% turn into 4x4 matrix
R_4mat = eye(4) ;
R_4mat(1:3, 1:3) = R ;

% convert matrix to tform
tform_r1 = affine3d(R_4mat) ;

% apply rotation (tform_init1) to segmented vnc
pc_seg_new = pctransform(pc_seg_new, tform_r1) ;

% store tform
tform_cell{end+1} = tform_r1 ;

% if debugFlag
%     figure ;
%     pcshowpair(pc_ref, pc_seg_new)
%     keyboard
% end
% -------------------------------------------------------------------------
% estimate scaling to get the two point clouds to be the right general size
seg_min = min(pc_seg_new.Location) ;
seg_max = max(pc_seg_new.Location) ;

ref_min = min(pc_ref.Location) ;
ref_max = max(pc_ref.Location) ;

s_val = (ref_max- ref_min)./(seg_max - seg_min) ;

% apply scaling
tform_s1 = affine3d([s_val(1), 0, 0, 0; 0, s_val(2), 0, 0 ; ...
    0, 0, s_val(3), 0 ; 0, 0, 0, 1]) ;
pc_seg_new = pctransform(pc_seg_new, tform_s1) ;

% store tform
tform_cell{end+1} = tform_s1 ;

% if debugFlag
%     figure ;
%     pcshowpair(pc_ref, pc_seg_new)
%     keyboard
% end

% ---------------------------------------------------------
% translate seg vnc to the cm of ref vnc
t_vec2 = mean(pc_ref.Location) - mean(pc_seg_new.Location) ;

% turn translation vector into tform
tform_t2 = affine3d([1, 0, 0, 0; 0, 1, 0, 0 ; 0, 0, 1, 0 ; ...
    t_vec2, 1]) ;

% apply translation
pc_seg_new = pctransform(pc_seg_new, tform_t2) ;

% store tform
tform_cell{end+1} = tform_t2 ;

if debugFlag
    figure ;
    pcshowpair(pc_ref, pc_seg_new)
    keyboard
end

% ------------------------------------------------------------------
%% calculate alpha shapes of point clouds to generate new images
% okay, this is getting convoluted, but what i think i'm finding is that
% the ref image has a lot more interior points. this leads to bad point
% cloud registration. since i only want the outside boundaries to match up,
% i'm going to get alpha shapes for each of these point clouds
fprintf('Finding alpha shapes for point clouds...\n')

% so first get alpha shapes
shp_ref = alphaShape(pc_ref.Location(:,1), pc_ref.Location(:,2),...
    pc_ref.Location(:,3)) ;
shp_seg = alphaShape(pc_seg_new.Location(:,1), pc_seg_new.Location(:,2),...
    pc_seg_new.Location(:,3)) ;

% then we're going to get the boundary facet points of these alpha shapes
[~, xyz_ref] = boundaryFacets(shp_ref) ;
[~, xyz_seg] = boundaryFacets(shp_seg) ;

if debugFlag
    figure ;
    hold on
    plot3(xyz_ref(:,1), xyz_ref(:,2), xyz_ref(:,3), 'b.')
    plot3(xyz_seg(:,1), xyz_seg(:,2), xyz_seg(:,3), 'r.')
    axis equal ; box on ; grid on
    keyboard
end

% make new point clouds from boundary facets
pc_ref_boundary = pointCloud(xyz_ref) ;
pc_seg_boundary = pointCloud(xyz_seg) ;
% ------------------------------------------------------------------
%% try point cloud registration
fprintf('Attempting registration...\n')
% down sample point clouds
fixed = pcdownsample(pc_ref_boundary, 'gridAverage',grid_ss_ref);
moving = pcdownsample(pc_seg_boundary, 'gridAverage', grid_ss_seg);

if debugFlag
    figure ;
    pcshowpair(fixed, moving)
    keyboard
end
% now perfrom PC registration
tform_reg = pcregistercpd(moving,fixed, 'Transform','rigid',...
    'MaxIterations', 200, 'Verbose',true);

% align seg point cloud and make plot
pc_aligned = pctransform(moving, tform_reg);

% plot results?
if plotFlag
    figure ;
    pcshowpair(fixed,pc_aligned)
end

% save tform?
tform_cell{end+1} = tform_reg ;
savePathFull = fullfile(savePath, 'tform_cell.mat') ;
if ~exist(savePathFull,'file') || overWriteFlag
    save(savePathFull,'tform_cell')
else
    fprintf('Overwrite old tform? \n')
    keyboard
end
% ------------------------------------------------------------------------
%
%{
% test transformation?
test = imwarp(segVNC, tform_init, 'OutputView', imref3d(size(refVNC_gray))) ;
% ---------------------------------------------------------
%% perform image registration
fprintf('Attempting registration...\n')
[optimizer,metric] = imregconfig('monomodal') ;
tform = imregtform(segVNC_gray,refVNC_gray,'affine',optimizer,metric, ...
    'DisplayOptimization', true, 'InitialTransformation', tform_init) ;

% test registration
segVNC_reg = imwarp(segVNC_gray,tform,...
    'OutputView',imref3d(size(refVNC_gray)));

figure ; imshowpair(max(refVNC_gray,[],3), max(segVNC_reg,[],3)) ;
%}
