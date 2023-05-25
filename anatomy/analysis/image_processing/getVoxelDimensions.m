% -------------------------------------------------------------------------
% quick function to grab voxel dimensions (i.e. vector [sx, sy, sz] where
% si is the size in microns of a voxel in the ith dimension
% -------------------------------------------------------------------------
function voxDimVec = getVoxelDimensions(dataFilename) 

% get image info
imgInfo = imfinfo(dataFilename) ; 

% keep this in try/catch because inconsistencies in image formating can
% make this difficult
try
    % x and y resolutions should be stored as fields of the imgInfo struct
    % for each image. NB: sometimes listed units are off, but i think we
    % can reasonably assume it's in microns
    voxSizeX = 1/(mean([imgInfo(:).XResolution])) ; 
    voxSizeY = 1/(mean([imgInfo(:).YResolution])) ; 
    
    % z resolution is given as a string in the ImageDescription field
    % (along with other potentially useful stuff)
    imgDescription = imgInfo(1).ImageDescription ; 
    search_str = 'spacing=' ; 
    imgDescriptionSplit = strsplit(imgDescription) ; 
    search_idx = cellfun(@(y) contains(y, search_str), imgDescriptionSplit) ; 
    
    if (sum(search_idx) == 1)
        spacing_str = imgDescriptionSplit{search_idx} ; 
        [~, end_idx] = regexp(spacing_str, search_str) ; 
        voxSizeZ = str2double(spacing_str((end_idx+1):end)) ;
    else
        % need to figure out what to do in this case
       voxSizeZ = nan ;  
    end
    voxDimVec = [voxSizeX, voxSizeY, voxSizeZ] ; 
catch
    voxDimVec = nan(1,3) ;
end


end