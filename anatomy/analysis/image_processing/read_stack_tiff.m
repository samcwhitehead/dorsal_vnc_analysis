%==========================================================================
% function to read in confocal stacks in tiff format
%
%   INPUTS:
%       =dataRoot = root folder where stacks are stored
%       -lineID = string identifying the fly line, e.g. 'BJD_SS00865'
%       -template_sex = string to determine whether to use 'FEMALE',
%       'MALE', or 'UNISEX' template
%       -CNS_region = string to determine whether we want 'VNC' or 'HR'
%       -file_type = file extention. should be tiff, but I guess it's good
%       to keep it flexible
%==========================================================================
function tiffMat = read_stack_tiff(dataRoot, lineID, template_sex, ...
    CNS_region, file_type) 
%------------------------------------------------
% deal with function inputs
if ~exist('template_sex','var')
    template_sex = 'UNISEX' ; 
end
if ~exist('CNS_region','var')
    CNS_region = 'VNC' ; 
end
if ~exist('file_type','var')
    file_type = '.tif' ; 
end
N_channels = 3 ; % can i read this automatically somehow?
%-----------------------------------------------
% information about data path
dataDir = dir(dataRoot) ; 
dataDir = dataDir(3:end) ; 

% get folder corresponding to line ID
folder_ind = arrayfun(@(x) contains(x.name, lineID),dataDir) ; 
if sum(folder_ind) ~= 1
    disp('could not find unique folder')
    keyboard 
end

dataFolder = fullfile(dataRoot,dataDir(folder_ind).name) ; 
dataFolderDir = dir(dataFolder) ; 

% get filename for appropriate image
filename_ind = arrayfun(@(x) contains(x.name,template_sex) & ...
    contains(x.name,CNS_region) & endsWith(x.name,file_type),...
    dataFolderDir) ;

if sum(filename_ind) ~= 1
    disp('could not find unique filename')
    keyboard
end

dataFilename = fullfile(dataFolder, ...
    dataFolderDir(filename_ind).name) ; 

%-----------------------------------------------
% get image info and initialize array to store stack
imgInfo = imfinfo(dataFilename) ; 

imHeight = imgInfo(1).Height ; 
imWidth = imgInfo(1).Width ; 
imBitDepth = imgInfo(1).BitDepth ;
N_images = length(imgInfo) ; 
imDepth = round(N_images/N_channels) ; 

switch imBitDepth
    case 16
        tiffMat = uint16(zeros(imHeight, imWidth, 3, imDepth)) ;
    case 8 
        tiffMat = uint8(zeros(imHeight, imWidth, 3, imDepth)) ;
    otherwise
        tiffMat = zeros(imHeight, imWidth, 3, imDepth) ;
end

%----------------------------------------------
% now load in images 
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