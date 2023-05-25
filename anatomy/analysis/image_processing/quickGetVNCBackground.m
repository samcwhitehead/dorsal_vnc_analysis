% -------------------------------------------------------------------------
% quick and dirty script to get VNC binary image as aggregate from
% old DN images
%
% should only need to do this once, but first need to generate our template
% VNC. going to sum all VNC bw images aligned to template and take average
% -------------------------------------------------------------------------
% get paths to data and where to save
dataPath = 'D:\Fly Imaging\Erica Interneuron Stacks\DN_old\vnc_templates\' ;
savePath = 'D:\Fly Imaging\Erica Interneuron Stacks\VNC_for_drawings\' ; 
dataDir = dir(fullfile(dataPath,'*_bw.mat')) ;

% overWriteFlag = true ; 

% shrinking image makes the boundary faster to calculate
shrinkFactor = 7 ; 

% loop over bw images
for k = 1:length(dataDir)
    % load current image
    vncBW = importdata(fullfile(dataDir(k).folder, dataDir(k).name)) ;
    
    % add image to sum
    if k == 1
        vncMat = double(vncBW) ;
    else
        vncMat = vncMat + double(vncBW) ;
    end
    
    disp(k)
end

% divide by number of images
vncMat = vncMat./length(dataDir) ;

% binarize vnc image
level = graythresh(vncMat) ; 
vncBW = imbinarize(vncMat, level) ;

% fill holes in image
vncBW = bwmorph3(vncBW, 'fill') ; 

% save vnc image and bw
save(fullfile(savePath, 'vncIm.mat'), 'vncMat') 
save(fullfile(savePath, 'vncBW.mat'), 'vncBW') 

% ------------------------------------------------------
% get boundary points of VNC for patch drawing
vncPerim = bwperim(vncBW) ;
vncPerim_smooth = imbinarize(smooth3(vncPerim, 'gaussian',9)) ;
nv = reducevolume(vncPerim_smooth, shrinkFactor) ;

[qy, qx, qz] = ind2sub(size(nv), find(nv(:))) ;
% image <-> cartesian
% qy = size(nv,1) - qy + 1 ;
%qz = size(nv,3) - qz + 1 ;
vnc_bound = boundary(qx, qy, qz) ;

% add boundary, points, and shrink factor to structure
vnc_struct = struct() ; 

vnc_struct.qx = qx ;
vnc_struct.qy = qy ; 
vnc_struct.qz = qz ;
vnc_struct.vnc_bound = vnc_bound ;

vnc_struct.shrinkFactor = shrinkFactor ; 
vnc_struct.vncBW = vncBW ; 

% save vnc boundary structure
save(fullfile(savePath, 'vnc_struct.mat'),'vnc_struct') 