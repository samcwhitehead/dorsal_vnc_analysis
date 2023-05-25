%--------------------------------------------------------------------------
% script to generate plot for comparing pulses (fast and slow) for
% experimental vs control flies
%--------------------------------------------------------------------------
%% PATH INFO AND PARAMS
% path to data file 
[mfilePath, ~, ~] = fileparts(mfilename('fullpath')) ; 
parentDirectory = fileparts(mfilePath) ;
dataPath = fullfile(parentDirectory, 'data') ;
savePath = fullfile(parentDirectory,'figs','output') ;

% define the lines/variables that we want to plot (currently limited by
% what's in data structure
MN_cell = {'DVM', 'tp2' ,  'ps1', 'i2' , 'hg1', 'hg2', 'hg3'} ;  % which motor neuron drivers to plot
driver_cell = {'SS41068', 'SS47120', 'SS47152', 'SS37246', 'SS48311',...
    'SS37253','SS49039'} ; % for the moment, pick out good examples of each MN
var_cell = {'pulseFrac','slow','fast'} ;
varName_cell = {'P_{slow} frac.', 'Slow\nPulse', 'Fast\nPulse'} ;

silencer_cell = {'Kir'} ;
exprCond_cell = {'ctrl','expr'} ;

% save output?
if ~exist(savePath,'dir')
    mkdir(savePath)
end
saveFlag = false ;

% general params
define_constants_kf ;

singFracMin = 5e-3 ;
singPercThresh = 5e-3 ;
minPulseNum = 20 ;

N_BOOT = 500 ;
Fs = 10000 ;
xlim_t = 10*[-1, 1] ;

N_plot_rows = length(var_cell) ;
N_plot_cols = length(driver_cell) ;

%--------------------------------------------------------------------------
%% LOAD DATA STRUCTURE
summaryStructPath = fullfile(dataPath, 'pulseSummaryStruct.mat') ;
pulseSummaryStruct = importdata(summaryStructPath) ; 

%--------------------------------------------------------------------------
%% PLOT PREFERENCES

%--------------------------------------------
%plot preferences
figPosition = [0.5972    0.8403    5.1    2.25] ;

% pulse drawings
lw_thick = 1.0 ;
CI_alpha = 0.4 ;

% data plots
showMM = 'median' ;

% xlabels = cell(1,2) ;
% xlim = [0.5, 2.5] ;
ylim_cell = {[0, 1], 0.3*[-1, 1],  0.3*[-1, 1]} ;
x_ticks = 1:(length(silencer_cell)*length(exprCond_cell)) ;
xlim = [0.5, length(x_ticks) + 0.5] ;

% subplot dimensions
gap = [0.05, 0.04] ; %[0.03, 0.05] ; %[0.1 0.10] ;
marg_h = [0.02, 0.1] ; %[0.08 0.08] ;
marg_w = [0.11 , 0.01] ;

%------------------------------------------------
%% MAKE FIGURE
h_main = figure('PaperPositionMode','auto','Units','inches','Position',...
    figPosition) ;

for ii = 1:length(driver_cell)
    for kk = 1:length(var_cell)
        %subPlotInd
        ax = subtightplot(N_plot_rows, N_plot_cols, ...
            (kk-1)*N_plot_cols + ii, gap, marg_h, marg_w) ; 
        hold on
        cc = 1 ;
        
        var_curr = var_cell{kk} ;
        for jj = 1:length(silencer_cell)
            silencer_curr = silencer_cell{jj} ;
            for mm = 1:length(exprCond_cell)
                exprCond_curr = exprCond_cell{mm} ;
                data_curr = pulseSummaryStruct(ii).([var_curr '_' silencer_curr ...
                    '_' exprCond_curr]) ;
                
                % color for given plot element
                if strcmp(exprCond_curr,'expr')
                    %                     cmap_curr = brewermap(6,cmap_struct.(driver_cell{ii})) ;
                    cmap_curr = brewermap(6,cmap_struct.SS37246) ;
                else
                    cmap_curr = brewermap(6,'Greys') ;
                end
                
                %----------------------------------------------------------
                %% if scatter plot...
                if size(data_curr,1) == 1
                    % get pVal
                    pVal = ranksum(pulseSummaryStruct(ii).([var_curr '_' silencer_curr ...
                        '_ctrl']), pulseSummaryStruct(ii).([var_curr '_' ...
                        silencer_curr '_expr'])) ;
                    
                    % adjust color by pvalue
                    if ((pVal <= 0.05) && (mm == 2)) || (mm == 1)
                        plotColor = cmap_curr(5,:)  ;
                    else
                        plotColor = [0, 0, 0] ;
                    end
                    
                    % plot data
                    ax = myPlotGroupDataToAxis(ax, cc, data_curr, ...
                        'showMM', showMM, ...
                        'plotColor', plotColor, ...
                        'barColor', plotColor) ;

                    % ------------------
                    % add sig stars?
                    % add sig star
                    if pVal < 0.001
                        text(0.5, 0.8, '***','Units','normalized')
                    elseif pVal < 0.01
                        text(0.5, 0.8, '**','Units','normalized')
                    elseif pVal < 0.05
                        text(0.5, 0.8, '*','Units','normalized')
                    end

                    %------------------------------------------------------
                    %% axis properties
                    box off
                    set(gca,'fontName',fontName,'fontSize',axisFontSize)
                    set(gca,'ylim',ylim_cell{kk})
                    % sort out y labels
                    if ii == 1
                        ylabel(varName_cell{kk}, 'fontName',fontName,...
                            'fontSize',labelFontSize)
                    else
                        set(ax,'YTickLabel',[])
                    end
                    
                   
                    set(ax,'xlim',xlim)
                    set(ax,'XTick',x_ticks,'XTickLabel',...
                        {'ctrl',silencer_curr})
                    ax = b1_paper_prettify_axis(ax, true, false, false, ...
                        false, true, 1.0 ) ;
                %----------------------------------------------------------
                %% if time series plot
                else
                    % get color for time series
                    plotColor = cmap_curr(end,:)  ;
                    
                    % get time info
                    pulse_t = size(data_curr,2) ;
                    T = 1000*((1:pulse_t) - round((pulse_t - 1)/2) - 1)/Fs ;
                    
                    % get data mean and CI
                    data_grand_mean = nanmean(data_curr) ;
                    CI = bootci(N_BOOT,@nanmean,data_curr) ;
                    
                    % plot CI
                    h_fill = fill(ax, [T,fliplr(T)],...
                        [CI(1,:),fliplr(CI(2,:))],...
                        plotColor,'linestyle','none');
                    h_fill.FaceAlpha = CI_alpha ;
                    % plot grand mean
                    plot(ax, T, data_grand_mean,'-','Color', plotColor,...
                        'LineWidth',lw_thick)
                    
                    set(gca,'xlim',xlim_t)
                    set(gca,'ylim',ylim_cell{kk})
                    
                    % remove axis lines
                    set(ax,'Color','none')
                    set(ax,'xcolor','none')
                    set(ax,'ycolor','none')
                    
                    % add y label if on the far left
                    if ii == 1
                        % make y label visible
                        ax.YAxis.Label.Color=[0 0 0]; 
                        ax.YAxis.Label.Visible='on';
                        
                        % set y label
                        ylab = ylabel(ax, sprintf(varName_cell{kk}), ...
                            'fontName',fontName,...
                            'fontSize',labelFontSize) ;
                        
                        % move y label closer to axis
                        set(ylab,'Units','normalized')
                        ylab_pos = get(ylab,'Position') ;
                        ylab_pos(1) = 0.5*ylab_pos(1) ;
                        set(ylab,'Position',ylab_pos)
                    else
                        set(ax,'YTickLabel',[])
                    end
                    
                    % ----------------------------------------
                    % add scale bar (to one plot)
                    if (ii == 1) && (kk == 2)
                        % draw line
                        sbx1 = -7.5 ;
                        sbx2 = -2.5 ;
                        sby1 = -0.2 ;
                        sby2 = sby1 ;
                        
                        plot(ax, [sbx1, sbx2], [sby1, sby2], 'k-', ...
                            'LineWidth', 1.5)
                        
                        % add text 
                        sb_txt = text(mean([sbx1, sbx2]), sby1 - 0.025, ...
                            sprintf('%d ms', sbx2 - sbx1), ...
                            'HorizontalAlignment', 'center',...
                            'VerticalAlignment', 'top',...
                            'fontName', fontName,...
                            'fontSize', labelFontSize) ; 
                    end
                    
                    % ----------------------------
                    % move 3rd row up
                    if kk == 3 
                        ax_pos = get(ax, 'Position') ;
                        ax_pos(2) = ax_pos(2) + 0.03 ;
                        set(ax, 'Position', ax_pos) 
                    end
                end
                cc = cc + 1 ;
            end
            
        end
        
    end
end

if saveFlag
    set(h_main, 'Renderer', 'painters')
    savefig(h_main,fullfile(savePath,'pulseSummary.fig'))
    print(h_main,fullfile(savePath,'pulseSummary.svg'),'-dsvg')
end
