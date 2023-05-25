%==========================================================================
% generate maximum intensity projection images for the 3 channels of a
% confocal stack
%==========================================================================
function h_mips = plot_channel_mips(imMat)
% generate single-color colormaps
N_colors = 500 ;
blue_cmap = [zeros(N_colors,2), linspace(0,1,N_colors)'] ;
red_cmap =[ linspace(0,1,N_colors)', zeros(N_colors,2)] ;
green_cmap = [zeros(N_colors,1), linspace(0,1,N_colors)', ...
    zeros(N_colors,1)] ;
gray_cmap = repmat(linspace(0,1,N_colors)',1,3) ;

% define some useful constants and figure out color order
defineStacksConstants

if size(size(imMat),2) == 4
    % this assumes full rgb tiff stack
    plot_order = [BLUE, RED, GREEN] ;
    plot_titles = {'Blue Channel', 'Red Channel', 'Green Channel'} ;
    cmap_cell = {blue_cmap, red_cmap, green_cmap} ;
elseif size(size(imMat),2) == 3
    % in case the input is just a single channel, account for that
    plot_order = 1 ;
    plot_titles = {'Channel'} ;
    cmap_cell = {hot} ;
else
    error('Channel dimensions are off')
end

% make figure
h_mips = figure ;
ax_arr = gobjects(length(plot_order),1) ;
for i = plot_order
    ax_arr(i) = subplot(1,length(plot_order),i) ;
    mip_curr = make_channel_mip(imMat,i) ;
    imshow(mip_curr,[])
    colormap(ax_arr(i), cmap_cell{i}) ;
    title(plot_titles{i})
    
end

end

%--------------------------------------------------------------------------
% accessory function to make mip
function mip_out = make_channel_mip(imMat, channel_num)

if size(size(imMat),2) == 4
    mip_out = squeeze(max(imMat(:,:,channel_num,:),[],4)) ;
elseif size(size(imMat),2) == 3
    mip_out = squeeze(max(imMat(:,:,:),[],3)) ;
else
    error('Channel dimensions are off')
end

end