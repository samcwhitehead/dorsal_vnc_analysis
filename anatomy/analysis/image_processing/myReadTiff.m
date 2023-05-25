% -------------------------------------------------------------------------
% function to (quickly?) read in tiff images as a 8bit or 16bit int matrix
% -------------------------------------------------------------------------
function [tiffMat, voxDimVec] = myReadTiff(dataFilename, N_channels)
% ----------------------------
% params and inputs
if ~exist('N_channels','var') || isempty(N_channels)
   N_channels = 3 ;  
end

% ------------------------------------------------
% get image info
imgInfo = imfinfo(dataFilename) ; 

imHeight = imgInfo(1).Height ; 
imWidth = imgInfo(1).Width ; 
imBitDepth = imgInfo(1).BitDepth ;
N_images = length(imgInfo) ; 
imDepth = round(N_images/N_channels) ; 

switch imBitDepth
    case 16
        tiffMat = uint16(zeros(imHeight, imWidth, N_channels, imDepth)) ;
    case 8 
        tiffMat = uint8(zeros(imHeight, imWidth, N_channels, imDepth)) ;
    otherwise
        tiffMat = zeros(imHeight, imWidth, N_channels, imDepth) ;
end

% including the voxel dimensions, in units of microns. if this fails, just
% returns nan values
voxDimVec = getVoxelDimensions(dataFilename)  ;
%----------------------------------------------
% now load in images 
switch N_channels
    case 1
        tiffMat = squeeze(tiffMat) ; 
        for i = 1:N_images
            tiffMat(:,:,i) = imread(dataFilename,i) ; 
        end
    otherwise
        TifLink = Tiff(dataFilename,'r') ;
        cc = 1 ;
        for i = 1:N_images
            channel_curr = mod(i,N_channels) + 1 ;
            TifLink.setDirectory(i) ;
            tiffMat(:,:,channel_curr,cc) = TifLink.read() ;
            if channel_curr == 2
                cc = cc + 1 ;
            end
        end
end

end