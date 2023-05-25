%--------------------------------------------------------------------------
% script to generate some summary stats plots for courtship song data
%
% change from v3: option for plotting median + iqr or mean +/- se
%--------------------------------------------------------------------------
%% PATH INFO AND PARAMS
% path to data file 
[mfilePath, ~, ~] = fileparts(mfilename('fullpath')) ; 
parentDirectory = fileparts(mfilePath) ;
dataPath = fullfile(parentDirectory, 'data') ;
savePath = fullfile(parentDirectory,'figs','output') ;

% define the lines/variables that we want to plot (currently limited by
% what's in data structure
MN_cell = {'DVM', 'tp2' ,  'ps1', 'i2' , 'hg1', 'hg2', 'hg3'} ; % , 'hg1 + hg3/4', 'b1 + hg4?', 'hg1 + i1', 'i1 + i2'} ; % which motor neuron drivers to plot
driver_cell = {'SS41068', 'SS47120', 'SS47152', 'SS37246', 'SS48311',...
    'SS37253','SS49039'} ; % , 'SS40864', 'SS40980', 'SS45772', 'SS45782'} ; % for the moment, pick out good examples of each MN
var_cell = {'fracSinging','pulseIndex','sineIndex'} ;
varName_cell = {'% Singing', 'Pulse Index', 'Sine Index'} ;

silencer_cell = {'Kir'} ;
exprCond_cell = {'ctrl','expr'} ;

define_constants_kf ;

% save results?
saveFlag = false ;
if ~exist(savePath,'dir')
    mkdir(savePath)
end

% params
singFracMin = 5e-3 ; %5e-3 ;
singPercThresh = singFracMin ; % 5e-3 ;
N_plot_rows = length(var_cell) ;
N_plot_cols = length(driver_cell) ;

% ------------------------------------------------------------------------
%% LOAD DATA
songDataStruct = importdata(fullfile(dataPath, 'songDataStruct.mat')) ; 

% ------------------------------------------------------------------------
%% SET PLOT PREFERENCES

figPosition = [0.5972    0.8403    5.1    3.8] ;

% bar plots
bar_width = 0.7 ; % 0.5
barEdgeColor = 'k' ; 
barAlphaVal = 0.33 ;

% data plots
showMM = 'median' ; 

cmap_grays = brewermap(6,'Greys') ; 
% xlabels = cell(1,2) ;
% xlim = [0.5, 2.5] ;
ylim_cell = {[0, 100], [0, 0.15], [0, 0.4]} ;
x_ticks = 1:(length(silencer_cell)*length(exprCond_cell)) ;
xlim = [0.5, length(x_ticks) + 0.5] ;

% subplot dimensions
gap = [0.06, 0.04] ; %[0.03, 0.05] ; %[0.1 0.10] ;
marg_h = [0.05, 0.1] ; %[0.08 0.08] ;
marg_w = [0.075, 0.01] ;

%------------------------------------------------
%% MAKE FIGURE
h_main = figure('PaperPositionMode','auto','Units','inches','OuterPosition',...
    figPosition) ;

% ---------------------------------------
% loop over drivers
for ii = 1:length(driver_cell)
    % ------------------------------------
    % loop over variables
    for kk = 1:length(var_cell)
        %subPlotInd
        ax = subtightplot(N_plot_rows, N_plot_cols, ...
            (kk-1)*N_plot_cols + ii, gap, marg_h, marg_w) ;
        hold on
        cc = 1 ;
        
        var_curr = var_cell{kk} ;
        % -------------------------------------------------
        % loop over effectors
        for jj = 1:length(silencer_cell)
            silencer_curr = silencer_cell{jj} ;
            
            % --------------------------------------------------------
            %% get p value
            if strcmp(var_curr,'fracSinging')
                [~,pVal,~] = fishertest(songDataStruct(ii).SongFracTable) ;
            else
                pVal = ranksum(songDataStruct(ii).([var_curr '_' silencer_curr ...
                    '_ctrl']), songDataStruct(ii).([var_curr '_' ...
                    silencer_curr '_expr'])) ;
            end
            
            % multiple hypothesis correction
            pVal = length(var_cell)*pVal ; 
            % ----------------------------------------------
            % loop over conditions
            for mm = 1:length(exprCond_cell)
                exprCond_curr = exprCond_cell{mm} ;
                data_curr = songDataStruct(ii).([var_curr '_' silencer_curr ...
                    '_' exprCond_curr]) ;
                
                % color for given plot element
                if strcmp(exprCond_curr,'expr')
                    cmap_curr = brewermap(6,cmap_struct.SS37246) ;
                else
                    cmap_curr = brewermap(6,'Greys') ;
                end
                
                % check if this data needs a bar or scatter plot
                if length(data_curr) == 1
                    
                    % change coloring based on whether or not it's
                    % significant
                    if ((pVal <= 0.05) && (mm == 2)) || (mm == 1)
                        plotColor = cmap_curr(end,:) ;
                    else 
                        plotColor = 'none' ;
                    end
                    
                    plt_curr = histogram(ax,...
                        'BinEdges',cc + (bar_width/2).*[-1,1],...
                        'BinCounts',data_curr,...
                        'FaceColor',plotColor,'EdgeColor',barEdgeColor,...
                        'FaceAlpha',barAlphaVal) ;
                    %set(get(plt_curr,'Children'),'FaceAlpha',barAlphaVal)
                else
                    
                    
                    % change coloring and filling based on whether or not it's
                    % significant
                    if ((pVal <= 0.05) && (mm == 2)) || (mm == 1)
                        plotColor = cmap_curr(5,:) ;
                    else 
                        plotColor = 'k' ;
                    end
                    
                    % --------------------------------------------
                    % plot data
                    ax = myPlotGroupDataToAxis(ax, mm, data_curr, ...
                        'showMM', showMM, ...
                        'plotColor', plotColor, ...
                        'barColor', plotColor) ;
                end
                cc = cc + 1 ;
            end
            
        end
        
        set(gca,'fontName',fontName,'fontSize',axisFontSize)
        set(gca,'xlim',xlim,'ylim',ylim_cell{kk})
        
        % sort out x labels
        if kk == length(var_cell)
            set(gca,'XTick',x_ticks,'XTickLabel',{'ctrl',silencer_curr})
        else
            set(gca,'XTick',x_ticks,'XTickLabel',[])
        end
        
        % sort out y labels
        if ii == 1
            ylabel(varName_cell{kk}, 'fontName',fontName,...
                'fontSize',labelFontSize)
        else
            set(gca,'YTickLabel',[])
        end
        
        % sort out title
        if kk == 1
            title(sprintf('%s\n({{\\bf%s}})',driver_cell{ii},MN_cell{ii}),...
                'fontName',fontName,'fontSize',titleFontSize,...
                'fontWeight','normal')
%             title(MN_cell{ii},'fontName',fontName,'fontSize',titleFontSize)
        end
        
        % add sig star
        if pVal < 0.001
            text(0.5, 0.8, '***','Units','normalized')
        elseif pVal < 0.01
            text(0.5, 0.8, '**','Units','normalized')
        elseif pVal < 0.05
            text(0.5, 0.8, '*','Units','normalized')
        end
        
        ax = gca ; 
        ax = prettify_axis(ax, true, false, false, false,...
            true, 1.0 ) ; 
        tit = get(ax,'title') ;
        % set(tit,'fontWeight','bold')
        hold off
    end
end

if saveFlag
    set(h_main, 'Renderer', 'painters')
    savefig(h_main,fullfile(savePath,'songIndexStats.fig'))
    print(h_main,fullfile(savePath,'songIndexStats.svg'),'-dsvg')
end

