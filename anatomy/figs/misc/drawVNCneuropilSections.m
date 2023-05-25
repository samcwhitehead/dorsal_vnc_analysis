% -------------------------------------------------------------------------
% function to generate drawing of full VNC with neuropil regions colored
% -------------------------------------------------------------------------
function h_main = drawVNCneuropilSections(regionLabels, colorVec, ...
    rootPath, regionAlpha, materialType)
%--------------------------
%% params and inputs
if ~exist('regionLabels','var') || isempty(regionLabels)
    regionLabels = {'NTct', 'WTct', 'HTct', 'IntTct','LTct', 'T1LegNp',...
        'T2LegNp', 'T3LegNp', 'AMNp', 'ANm', 'VAC1', 'VAC2', 'VAC3',...
        'mVAC'} ;
end
if ~exist('colorVec','var') || isempty(colorVec)
    colorVec = brewermap(length(regionLabels), 'Set1') ; %lines(length(regionLabels)) ; %0.6*ones(length(regionLabels),3) ; 
end
if ~exist('rootPath','var') || isempty(rootPath)
    rootPath = ['D:\Fly Imaging\Erica Interneuron Stacks\VNC neuropil\'...
        'segmented_images\registered\'] ; 
end
if ~exist('regionAlpha','var') || isempty(regionAlpha)
    regionAlpha = 1.0 ; %0.9 ; 
end
if ~exist('materialType','var') || isempty(materialType)
    materialType = 'metal' ; % 'dull' | 'shiny' 
end
scale = [1.0, 1.0, 1.0] ;
coord_z_scale = 1.0 ; %1.5 ; 

VNC_color = 0.9*[1 1 1] ; 
VNC_alpha = 0.5 ; 

az = 180 ; 
el = 90 ; 
% -------------------------
%% initialize figure 
h_main = figure ; 
hold on
ax = gca ; 
parent = hgtransform('Parent',ax);

% ---------------------------------------
%% first draw VNC
vnc_data = load(fullfile(rootPath, 'VNC.mat')) ; 
% want to scale VNC data to be slightly larger than other regions, thus
% avoiding overlap
vnc_coords = vnc_data.vox_coord_cell{1} ; 
vnc_coords_centroid = nanmean(vnc_coords) ; 
vnc_coords = vnc_coords - vnc_coords_centroid ;
vnc_coords(:,3) = vnc_coords(:,3) + 2 ; %1.2* + 3
vnc_coords = 1.05.*vnc_coords + vnc_coords_centroid ; 
VNC_grp = draw3Dboundary(ax, parent, vnc_coords, ...
    vnc_data.vox_boundary_cell{1}, VNC_color, VNC_alpha, scale) ;
% --------------------------------------
%% loop through regions
for i = 1:length(regionLabels)
   labelCurr = regionLabels{i} ; 
   % load current data
   vox_data = load(fullfile(rootPath, [labelCurr '.mat'])) ; 
   
   numObjects = length(vox_data.vox_coord_cell) ; 
   
   if strcmp(labelCurr,'VAC1') 
      translationVec = [0, 0, -10] ; 
   elseif strcmp(labelCurr,'VAC2') || strcmp(labelCurr,'VAC3')
       translationVec = [0, 0, -7] ; 
   else
       translationVec = [0, 0, 0] ; 
   end
   % ---------------------------------
   % loop through objects
   for j = 1:numObjects
       % scale z direction (doing it via hgtform seems to cause issues since
       % the data is off-center)
       vox_coords = vox_data.vox_coord_cell{j} ;
       vox_coords_mean = nanmean(vox_coords) ;
       vox_coords = vox_coords - vox_coords_mean ;
       vox_coords(:,3) = coord_z_scale.*vox_coords(:,3) ;
       vox_coords = vox_coords + vox_coords_mean ;
       
       % draw surface
       grp = draw3Dboundary(ax, parent, vox_coords, ...
           vox_data.vox_boundary_cell{j}, colorVec(i,:), regionAlpha, ...
           scale, translationVec, materialType) ;
   end
end

% ------------------------------------
%% set axis properties
axis equal 
axis tight
hlight = light ;
axis off
view(az, el)
end