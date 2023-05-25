% -------------------------------------------------------------------------
% function to make a plot of anatomical overlap of interneurons
% -------------------------------------------------------------------------
function [h_overlap, ax] = plot_overlap_data(overlap_data, labels_x, ...
    labels_y, plotStyle, normFlag, cmap)
% -------------------------
%% params and inputs
if ~exist('labels_x','var') || isempty(labels_x)
    labels_x = [] ;
end
if ~exist('labels_y','var') || isempty(labels_y)
    labels_y = [] ;
end
if ~exist('plotStyle','var') || isempty(plotStyle)
    plotStyle = 'sym' ; % can be 'sym' (symmetric) or 'tri' (lower triangular)
end
if ~exist('normFlag','var') || isempty(normFlag)
    normFlag = false ; % normalize overlap mat to max?
end
% if ~exist('clusterFlag','var') || isempty(clusterFlag)
%     clusterFlag = false ; % sort overlap mat based on hierarchical clustering?
% end
if ~exist('cmap','var') || isempty(cmap)
    cmap = brewermap([],'Blues') ;
end
% -------------------
% plot preferences
%cmap_name = 'Blues' ; %'YlGnBu' ;
ax_color = 'k' ; %'k' ; %'w' ; % 'k' ;
%c_axlim = [0, 0.4] ; % color bar limits
axisFontSize = 6 ; % needs to be tiny for labels
axisFontSizeSmall = 4.5 ; % needs to be tiny for labels

[N_stacks1, N_stacks2] = size(overlap_data) ;

% axis position (dendrogram and matrix):
gap = [0.00, 0.00] ; %[0.1 0.10] ;
marg_h = [0.01, 0.10] ; %[0.1 0.02] ;
marg_w =  [0.15, 0.125] ; % [0.125, 0.125] ;

rmZerosFlag = true ;
% -----------------------------
%% process overlap matrix to accommodate plot type
% create either symmetric or lower tri matrix, depending on plotStyle
if strcmp(plotStyle,'tri')
    % fill out other lower triangular matrix
    for ii = 1:N_stacks1
        for jj = 1:ii
            overlap_data(ii,jj) = overlap_data(jj,ii) ;
            overlap_data(jj,ii) = nan ;
        end
    end
end

% ---------------------------------------------------------------
%% normalize overlap matrix to max val ?
% NB: in most of my code, normalization happens prior to this
if normFlag
    overlap_data = normalize_overlap_mat(overlap_data) ;
end

% ---------------------------------------------------------------
%% get overlap values (different for table vs mat)
if istable(overlap_data)
    overlap_vals = overlap_data{:,:} ;
else
    overlap_vals = overlap_data ;
end

% -----------------------------------------------------------
%% change zero values to nan so they get different color?
if rmZerosFlag
    % repace zeros with nan
%     zero_idx = (overlap_vals == 0) ;
%     overlap_vals(zero_idx) = nan ;
    
    % set axis color to white so it stands out
    ax_color = 'w' ; %k
end
% -----------------------------------------------------------
%% make figure
h_overlap = figure('PaperPositionMode','auto') ;
ax = subtightplot(1,1,1, gap, marg_h, marg_w) ;
hold on
% overlap_vals = log(overlap_vals) ;
imagesc(ax, flipud(overlap_vals), 'AlphaData',~isnan(flipud(overlap_vals)))
set(ax,'color',ax_color)
%set(ax,'Ydir','normal')
% cmap = [0.5, 0.5, 0.5 ; cmap] ;
colormap([[0,0,0]; cmap])

% axis limits
%axis equal
axis tight
%caxis(c_axlim)

% axis labels
set(ax,'fontName','arial','fontSize',axisFontSizeSmall,...
    'TickLabelInterpreter','none')
set(ax,'XTick',1:N_stacks2, 'XTickLabel',labels_x,'XTickLabelRotation',90)
set(ax,'YTick',1:N_stacks1, 'YTickLabel',flipud(labels_y))

% if (N_stacks2 > 100)
%     ax.XAxis.FontSize = axisFontSizeSmall ;
% end
% if (N_stacks1 > 100)
%     ax.YAxis.FontSize = axisFontSizeSmall ;
% end

% color bar
c = colorbar('Location','southoutside') ;
c.Label.String = 'normalized overlap' ;
% if strcmpi(normType, 'none')
%     c.Label.String = 'volume overlap (vox num)' ;
% else
    
% end

% remove tick lengths
set(ax,'TickLength',[0 0])


end