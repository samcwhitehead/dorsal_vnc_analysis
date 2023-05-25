% -------------------------------------------------------------------------
% quick script to calculate alpha shapes and bw images (masks) for
% segmented neuropil region data
% -------------------------------------------------------------------------
%% path info and params
dataRoot = 'D:\Fly Imaging\Erica Interneuron Stacks\VNC neuropil\' ;
registeredFlag = true ; % use registered images?
if registeredFlag
    dataPath = fullfile(dataRoot, 'segmented_images', 'registered') ;
else
    dataPath = fullfile(dataRoot, 'segmented_images') ;
end
savePath = dataPath ;
segDataPathFull = fullfile(dataRoot, 'seg.tif') ;

ignoreFnList = {'tform_cell.mat'} ; % any files that aren't neuropil data

% check that our coordinate conversions are okay?
debugFlag = false ;

% get directory of segemented neuropil data
neuropilDir = dir(fullfile(dataPath, '*.mat')) ;

% size for bw images
if registeredFlag
    imSize = [1119, 573, 333]; % this the size for all our binary images of neurons
else
    % if we're not using registered image, load original segmentation
    neuropilTiff = myReadTiff(segDataPathFull,1) ;
    imSize = size(neuropilTiff) ;
    
    % define pixel val to neuropil region correspondence
    pix2Label = {'LTct',... % 1
        'WTct', ... %2
        'AMNp', ... %3
        'mVAC', ... %4
        'VAC1', ... %5
        'NTct', ... %6
        'T1LegNp', ... %7
        'T2LegNp', ... %8
        'T3LegNp', ... %9
        'ANm', ... %10
        'HTct', ... %11
        'VAC2', ... %12
        'VAC3', ... %13
        'IntTct'  } ; %14
    
    % get unique pixel values
    pixVals = unique(neuropilTiff(:)) ;
    pixVals = pixVals(pixVals ~= 0) ; % exclude 0 (background)
end

% ------------------------------------------------------------
%% loop over neuropil files
for ind = 1:length(neuropilDir)
    %% read out path and filename
    % (reading out path is redundant but w/e)
    path = neuropilDir(ind).folder ;
    fn = neuropilDir(ind).name ;
    [~, np_name, ~] = fileparts(fn) ;
    
    fprintf('Processing %s ... \n', np_name)
    
    % check that we're looking at neuropil data
    if sum(ismember(fn, ignoreFnList)) > 0
        fprintf('Skipping %s \n', np_name)
        continue
    end
    
    % ---------------------------------------
    %% load/read neuropil data
    np_struct = importdata(fullfile(path, fn)) ;
    
    % read out data
    vox_coord_cell = np_struct.vox_coord_cell ;
    vox_boundary_cell = np_struct.vox_boundary_cell ;
    % we should already have binary images stored, but check just in case
    if isfield(np_struct, 'bw_cell')
        bw_cell = np_struct.bw_cell ;
        calcBWFlag = false ;
    else
        bw_cell = cell(length(vox_coord_cell),1) ;
        calcBWFlag = true ;
    end
    
    % find number of connected components in this neuropil image dataset
    N_cc = length(vox_coord_cell) ;
    
    % initialize storage for output
    shp_cell = cell(N_cc,1) ;
    
    % ---------------------------------------
    %% loop over objects
    for k = 1:N_cc
        % get current coords
        coords = vox_coord_cell{k} ;
        
        % calculate alpha shape based on coordinates
        shp = alphaShape(coords(:,1), coords(:,2), coords(:,3)) ;
        
        % add alpha shapeto output
        shp_cell{k} = shp ;
        
        % --------------------------------------------------------------
        %% get a binary image of neuopil?
        if calcBWFlag
            % if we're using registered images, going to try this by using
            % alpha shape
            if registeredFlag
                % get xyz limits for current coordinates; use this to create mesh
                min_vals = min(coords) ;
                max_vals = max(coords) ;
                
                [xx, yy, zz] = meshgrid(min_vals(1):max_vals(1), ...
                    min_vals(2):max_vals(2), ...
                    min_vals(3):max_vals(3)) ;
                
                mesh_pts = [xx(:), yy(:), zz(:)] ;
                clear xx yy zz
                
                % test which points in meshgrid are in alpha shape
                in_np_idx = inShape(shp, mesh_pts(:,1), mesh_pts(:,2), ...
                    mesh_pts(:,3));
                coords_in_shp = uint16(mesh_pts(in_np_idx,:)) ;
                
                % initialize storage
                bw = false(imSize) ;
                
                % get convert coordinates (sub) to image index
                bw_idx = sub2ind(imSize, coords_in_shp(:,2), ...
                    coords_in_shp(:,1), coords_in_shp(:,3)) ;
                
                % mark these voxels as "true"
                bw(bw_idx) = true ;
                
            else
                % otherwise use segmented image
                % get pixel value corresponding to current neuropil region
                pixValCurr = find(cellfun(@(y) strcmpi(np_name, y),...
                    pix2Label));
                if numel(pixValCurr) ~= 1
                    fprintf('Error: no match for neuropil name: %s\n',...
                        np_name)
                    continue
                end
                
                % generate bw image
                bw = false(imSize) ;
                bw_idx = find(neuropilTiff == pixValCurr) ;
                bw(bw_idx) = true ;
            end
            
            % ----------------------------------------
            %% check output?
            if debugFlag
                figure ;
                imshow(max(bw,[],3))
                hold on
                plot(coords(:,1), coords(:,2), 'r.')
                
                keyboard
            end

            bw_cell{k} = bw ;
        end
    end
    
    % -----------------------------------------------------
    %% save np data with new bits
    savePathFull = fullfile(savePath, fn) ;
    save(savePathFull, 'vox_coord_cell', 'vox_boundary_cell', ...
        'shp_cell', 'bw_cell')
    
    fprintf('Finished processing %s \n', np_name)
end