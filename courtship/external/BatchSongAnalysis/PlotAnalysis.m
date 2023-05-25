function PlotAnalysis(folders,genotypes,control_folders,control_genotypes,savePath)
%Plot all results from multiple genotypes both with and without outliers, 

if nargin < 5
    PlotResults(folders,genotypes,control_folders,control_genotypes,1);
    %PlotResults(folders,genotypes,control_folders,control_genotypes,0)
else
    PlotResults(folders,genotypes,control_folders,control_genotypes,1,savePath);
    %PlotResults(folders,genotypes,control_folders,control_genotypes,0,savePath)
end





function PlotResults(folders,genotypes,control_folders,control_genotypes,remove_outliers,savePath)

% folders = Cell array of folders containing genotypes to plot
% genotypes = Cell array of genotypes to plot
% control_folders = Cell array of folders containing control
% control_genotypes = Cell array of control genotype(s). All control
% genotypes are included in all plots, but statistics are calculated for
% grand mean of controls.
% remove_outliers = logical (1/0) to plot with outliers removed or original
%
% USAGE
% PlotAnalysis({'folder1' 'folder2'},{'one' 'two' 'three'},{'folder3' 'folder1'},{'Ore-R' 'Canton-S'})
%

%collect control data

if nargin < 5
    remove_outliers = 1;
    savePath = folders{1} ; 
end

numControls = numel(control_genotypes);
numGenotypes = numel(genotypes);
SampleSize = zeros(1,numControls+numGenotypes);


controls = struct();

for i  = 1:numControls
    fprintf('%s\n',control_genotypes{i});
    count = 0;
    for ii = 1:numel(control_folders)
        if control_folders{ii}(end) ~= '\'
            control_folders{ii} = [control_folders{ii} '\'];
        end
        dir_list = dir(control_folders{ii});
%        for g = 1:numControls %check on whether controls can be found
%            idx = ~cellfun('isempty',strfind({dir_list.name},control_genotypes{g}));
%            if g==1 && sum(idx) < 1
%                fprintf('Control genotype not found. Check details.\n')
%                return
%            end
 %       end
        
        for iii = 1:numel(dir_list)
            file = dir_list(iii).name;
            [~,~,ext] = fileparts(file);
            if ~isempty(strfind(file,control_genotypes{i})) && strcmp(ext,'.mat')
                count = count + 1;
                AR = load([control_folders{ii} file], '-mat');
                if count == 1
                    ControlResults = AR.Stats2Plot;
                else
                    ControlResults = [ControlResults AR.Stats2Plot];
                end
            end
        end
    end
    varname = genvarname(control_genotypes{i});
    controls.(varname) = ControlResults;
end

    
for x = 1:numControls
    varname = genvarname(control_genotypes{x});
    SampleSize(x) = numel(controls.(varname));
end

%collect target results

results = struct();
for i  = 1:numel(genotypes)
    fprintf('%s\n',genotypes{i});
    count = 0;
    for ii = 1:numel(folders)
        if folders{ii}(end) ~= '\'
            folders{ii} = [folders{ii} '\'];
        end
        dir_list = dir(folders{ii});
%         for g = 1:numGenotypes
%             idx = ~cellfun('isempty',strfind({dir_list.name},genotypes{g}));
%             if sum(idx) < 1
%                 fprintf('Genotype not found. Check details.\n')
%                 return
%             end
%         end
        for iii = 1:numel(dir_list)
            file = dir_list(iii).name;
            [~,~,ext] = fileparts(file);
            if ~isempty(strfind(file,genotypes{i})) && strcmp(ext,'.mat')
                count = count + 1;
                AR = load([folders{ii} file], '-mat');
                if count == 1
                    Results = AR.Stats2Plot;
                else
                    Results = [Results AR.Stats2Plot];
                end
            end
        end
    end
    varname = genvarname(genotypes{i});
    results.(varname) = Results;
    
end
for y= 1:numGenotypes
    varname=genvarname(genotypes{y});
    SampleSize(y+numControls) = numel(results.(varname));
end
numSamples = sum(SampleSize(end-numGenotypes+1:end));
maxSampleSize = max(SampleSize);
names = fieldnames(results.(varname));

% densitycolors = {'blue','red','green','magenta','cyan'};

% geno_varname = genvarname(genotypes{i});
% numSamples = numel(results.(geno_varname));
% SampleSize(end) = numSamples;
% maxSampleSize = max(SampleSize);
% names = fieldnames(results.(geno_varname));

    

%make arrays for plotting
%for i = 1:numel(genotypes) %for each genotype

clf;
ha = tight_subplot(8,4,.05,.06,[.05 .03]);
for j = 1:27 %all results except models
    Results2Plot = NaN(maxSampleSize,numControls+numGenotypes);
    Trait = names{j};
    numOutliers= [];
    dataOutliers = [];
    color = 'k';
    if ~strcmp(Trait,'Sine2PulseNorm') && ~strcmp(Trait,'ipiDist') && ~strcmp(Trait,'lombStats')
        
        %collect control data
        for k = 1:numControls
            control_varname = genvarname(control_genotypes{k});
            for m = 1:numel(controls.(control_varname))
                Results2Plot(m,k) = controls.(control_varname)(m).(Trait)(1);
            end
        end
        %collect results
        for k = 1:numel(genotypes)
            geno_varname = genvarname(genotypes{k});
            for n = 1:numel(results.(geno_varname))
                Results2Plot(n,numControls+k) = results.(geno_varname)(n).(Trait)(1);
            end
        end
        OutliersRemovedResults2Plot = NaN(size(Results2Plot));
        if remove_outliers == 1
            for q = 1:size(Results2Plot,2) %for each column
                if sum(~isnan(Results2Plot(:,q)))>0 %if some data in columns
                    [OutliersRemovedResults2Plot(:,q),controlidx,~] = deleteoutliers(Results2Plot(:,q),.05,1);
%                     numOutliers = numel(controlidx) - sum(isnan(Results2Plot(:,q)));
%                     if numOutliers == 0
%                         numOutliers = [];
%                     end
%                 else
%                     OutliersRemovedResults2Plot(:,q) = Results2Plot(:,q);
%                     numOutliers = [];
                end
            end
        end
%             if sum(~isnan(Results2Plot(:,end)))>0 %sum(isnan(Results2Plot(:,end))) ~= numel(Results2Plot(:,end))
%                 [OutliersRemovedResults2Plot(:,end),dataidx,~] = deleteoutliers(Results2Plot(:,end),.05,1);
%                 dataOutliers= numel(dataidx) - sum(sum(isnan(Results2Plot(:,end))));
%                 if dataOutliers == 0
%                     dataOutliers = [];
%                 end
%             else
%                 OutliersRemovedResults2Plot(:,end) = Results2Plot(:,end);
%                 dataOutliers = [];
%             end
%         else
%             OutliersRemovedResults2Plot = Results2Plot;
%             controlOutliers= [];
%             dataOutliers = [];
%         end
        
        
        
        colors = cell(numControls + numGenotypes,1);
        colors(1:numControls) = {'k'};
        %colors{end} = color;
        %determine whether results are sign diff from controls
        for i = numControls + 1:numGenotypes + numControls
            h = ttest2(reshape(OutliersRemovedResults2Plot(:,1:numControls),1,numel(OutliersRemovedResults2Plot(:,1:numControls)))',...
                OutliersRemovedResults2Plot(:,i),0.01);
            if h == 1
                %change color of results
                color = 'r';
            else
                color = 'b';
            end
            colors{i} = color;
        end
        
        
        
        %plot in new panel
        axes(ha(j))
        title(Trait)
        errorbarjitter(OutliersRemovedResults2Plot,ha(j),'Ave_Marker','+','Colors',colors,'color_opts',[1 1 1])
        set(gca,'XTick',[],'YTickLabelMode','auto','XColor',get(gca,'Color'))
        %if numel(numOutliers) > 0
        %text(0.15,-.5,'#Outliers=' ,'Units','normalized', 'interpreter', 'none')
        %end
%         if numel(dataOutliers) > 0
%             text(0.5,-.1,['#Outliers=' num2str(dataOutliers)],'Units','normalized', 'interpreter', 'none')
%         end
        if min(min(OutliersRemovedResults2Plot)) < max(max(OutliersRemovedResults2Plot));
            ylim(ha(j),[min(min(OutliersRemovedResults2Plot)) max(max(OutliersRemovedResults2Plot))]);
        end
        
    elseif strcmp(Trait,'Sine2PulseNorm')%plot scatter plots for normalized sine to pulse
        NormS2P2Plot = NaN(maxSampleSize,2*(numControls+numGenotypes));
        %collect control data
        for k = 1:numControls
            control_varname = genvarname(control_genotypes{k});
            for m = 1:numel(controls.(control_varname))
                p = k + (k-1);
                NormS2P2Plot(m,[p p+1]) = controls.(control_varname)(m).(Trait)([1 2]);%place 2 columns of data in neighboring columns
            end
        end
        %collect results
        
        for m = 1:numel(genotypes)
            geno_varname = genvarname(genotypes{m});
            for n = 1:numel(results.(geno_varname))
                p = k + m + (k+m-1);
                NormS2P2Plot(n,[p p+1]) = results.(geno_varname)(n).(Trait)([1 2]);%place 2 columns of data in neighboring columns
            end
        end

        OutliersRemovedNormS2P2Plot = NaN(size(NormS2P2Plot));
        if remove_outliers == 1
            for m = 1:size(NormS2P2Plot,2)
                if sum(~isnan(NormS2P2Plot(:,m)))>0 %sum(isnan(NormS2P2Plot(:,m))) ~= numel(NormS2P2Plot(:,m))
                    OutliersRemovedNormS2P2Plot(:,m) = deleteoutliers(NormS2P2Plot(:,m),.05,1);
                else
                    OutliersRemovedNormS2P2Plot(:,m) = NormS2P2Plot(:,m);
                    controlidx = [];
                end
            end
        end
        
        OutliersRemovedNormS2P2Plot = NormS2P2Plot;
        controlidx= [];
        dataidx = [];
        
        
        %determine whether results are sign diff from controls
        %make arrays for aoctool
        x = reshape(OutliersRemovedNormS2P2Plot(:,1:2:end-1),[],1);
        y = reshape(OutliersRemovedNormS2P2Plot(:,2:2:end),[],1);
        num = size(OutliersRemovedNormS2P2Plot,1);
        group=zeros(numel(x),1);
        group(num+1:end)=1;
        group = group(isfinite(x));
        x = x(isfinite(x));	%new code to keep only non-infinite data
        y = y(isfinite(y));
        %group = zeros(size(x,1),1);
        %group(end-numSamples:end) = 1;
        [~,T,h] = aoctool(x,y,group,[],[],[],[],'off');
        
        [GroupRow,~] = find(strcmp(T,'group') == 1);
        [GroupXRow,~] = find(strcmp(T,'group*x') == 1);
        [~,Probcol] = find(strcmp(T,'Prob>F') == 1);
        
        
        if T{GroupRow,Probcol} < 0.01 || T{GroupXRow,Probcol} < 0.01
            %change color of results
            color = 'r';
        else
            color = 'b';
        end
        
        %plot in new panel
        axes(ha(j))
        title(Trait)
        hold on
        for k = 1:numControls
            p = k + (k-1);
            x = OutliersRemovedNormS2P2Plot(:,p);
            y = OutliersRemovedNormS2P2Plot(:,p+1);
            scatter(x,y,'k')
            try
                brob = robustfit(x,y);
                plot(x,brob(1)+brob(2)*x,'k')
            catch
            end
        end
        for m = 1:numGenotypes
            p = k + m + (k+m-1);
            x = OutliersRemovedNormS2P2Plot(:,p);
            y = OutliersRemovedNormS2P2Plot(:,p+1);
            scatter(x,y,color)
            try
                brob = robustfit(x,y);
                plot(x,brob(1)+brob(2)*x,'k')
            catch
            end
        end
        hold off
        set(gca,'YTickLabelMode','auto')
        set(gca,'XTickLabelMode','auto')
        ylim(ha(j),[min(min(OutliersRemovedNormS2P2Plot)) max(max(OutliersRemovedNormS2P2Plot))]);
        color = 'k';
        controlidx= [];
        dataidx = [];
        
        
        
    elseif strcmp(Trait,'ipiDist')%plot kernel density estimates for IPIs
        %ipiDist2Plot = NaN(maxSampleSize,2*(numControls+1));
        %collect control data
        
        axes(ha(j))
        title(Trait)
        hold on
        
        %for each control genotype
        for k = 1:numControls
            
            %for each sample in control genotype
            control_varname = genvarname(control_genotypes{k});
            
            %collect ipiDist data and plot in black
            for m = 1:numel(controls.(control_varname))
                plot(controls.(control_varname)(m).ipiDist.xi,controls.(control_varname)(m).ipiDist.f,'Color',[.5 .5 .5])
                alpha(.5)
            end
        end
        %axis tight
        cmap = hsv(numGenotypes);
        %for each sample in genotype
        for k = 1:numGenotypes
            geno_varname= genvarname(genotypes{k});
            for m = 1:numel(results.(geno_varname))
                %collect ipiDist data and plot in new color
                plot(results.(geno_varname)(m).ipiDist.xi,results.(geno_varname)(m).ipiDist.f,'Color',cmap(k,:))
                alpha(.5)
            end
        end
        xlim([200 700])
        set(gca,'XTickLabel',{'20','30','40','50','60','70'})
        
    elseif strcmp(Trait,'lombStats')
        %collect control data
        lombStats2PlotControlsF = {};
        lombStats2PlotControlsAlpha = {};
        count = 0;
        for k = 1:numControls
            control_varname = genvarname(control_genotypes{k});
            for m = 1:numel(controls.(control_varname))
                count = count+1;
                lombStats2PlotControlsF{count} = controls.(control_varname)(m).(Trait).F;
                lombStats2PlotControlsAlpha{count} = controls.(control_varname)(m).(Trait).Alpha;
            end
        end
        lombStats2PlotControlsF = cell2mat(lombStats2PlotControlsF');
        lombStats2PlotControlsAlpha = cell2mat(lombStats2PlotControlsAlpha');
        
        %collect results
        lombStats2PlotResultsF = {};
        lombStats2PlotResultsAlpha = {};
        count = 0;
        for k = 1:numGenotypes
            geno_varname= genvarname(genotypes{k});
            for n = 1:numel(results.(geno_varname))
                count = count+1;
                lombStats2PlotResultsF{count} = results.(geno_varname)(n).(Trait).F;
                lombStats2PlotResultsAlpha{count} = results.(geno_varname)(n).(Trait).Alpha;
            end
        end
        lombStats2PlotResultsF = cell2mat(lombStats2PlotResultsF');
        lombStats2PlotResultsAlpha = cell2mat(lombStats2PlotResultsAlpha');
        
        axes(ha(j))
        title(Trait)
        hold on
        scatter(lombStats2PlotControlsF,lombStats2PlotControlsAlpha,'k')
        scatter(lombStats2PlotResultsF,lombStats2PlotResultsAlpha,'b')
        hold off
        set(gca,'YTickLabelMode','auto')
        set(gca,'XTickLabelMode','auto')
        set(gca,'YDir','reverse')
        set(gca,'YSCale','log')
        xlim(gca,[0 0.1])
        if min([lombStats2PlotControlsAlpha;lombStats2PlotResultsAlpha]) < max([lombStats2PlotControlsAlpha;lombStats2PlotResultsAlpha])
            ylim(ha(j),[min([lombStats2PlotControlsAlpha;lombStats2PlotResultsAlpha]) max([lombStats2PlotControlsAlpha;lombStats2PlotResultsAlpha])]);
        end
    end
    
end

%collect and align models
numModels = numSamples + 1;
t = cell(numModels,1);
y=0;
for i = 1:numGenotypes
    geno_varname= genvarname(genotypes{i});
    for n = 1:numel(results.(geno_varname))
        y = y+1;
        t{y} = results.(geno_varname)(n).PulseModels.NewMean;
    end
end

t{y+1} = results.(geno_varname)(n).PulseModels.OldMean;

%remove empty models

t = t(~cellfun('isempty',t));

max_length = max(cellfun(@length,t));
total_length = 2* max_length;
Z = zeros(size(t,1),total_length);
if numGenotypes >0
    for n=1:size(t,1)
        if ~isempty(t{n})
            X = t{n};
            T = length(X);
            [~,C] = max(abs(X));%get position of max power
            %flip model is strongest power is negative
            if X(C) <0
                X = -X;
            end
            %center on max power
            left_pad = max_length - C;  %i.e. ((total_length/2) - C)
            right_pad = total_length - T - left_pad;
            Z(n,:) = [zeros(1,left_pad) X zeros(1,right_pad)];
        end
    end
end
[Z,~] = alignpulses(Z,2);
%trim down models
first = find(sum(Z,1),1,'first');
last = find(sum(Z,1),1,'last');
Z = Z(:,first:last);
axes(ha(j+1))
plot(ha(j+1),Z(end,:)','Color',[.6 .6 .6],'LineWidth',8)
hold on
m=0;
for k = 1:numGenotypes
    geno_varname= genvarname(genotypes{k});
    for n = 1:numel(results.(geno_varname))
        m = m + 1;
        plot(ha(j+1),Z(m,:)','Color',cmap(k,:));
    end
end
hold off
axis(ha(j+1),'tight');
axis([ha(j+1) ha(j+2) ha(j+3)],'off');


%print useful information in final panel
axes(ha(j+2))

%collect controls
%conGenos = [];
% for a= 1:numControls
%     conGenos = [conGenos char(control_genotypes{a})];
% end
text(0,1,['Controls = ' control_genotypes], 'interpreter', 'none')

%collect control folders
% conFolders = [];
% for a= 1:numel(control_folders)
%     folder = regexp(control_folders{a},'/','split');
%     conFolders= [conFolders folder(end-1)];
% end
text(.25,1,['Control Folders = ' control_folders], 'interpreter', 'none')

axes(ha(j+3))
text(0,1,['Genotypes = ' genotypes], 'interpreter', 'none')
%collect data folders
% resFolders = [];
% for a= 1:numel(folders)
%     folder = regexp(folders{a},'/','split');
%     resFolders= [resFolders folder(end-1)];
% end
text(.25,1,['Analysis Folders = ' folders], 'interpreter', 'none')

%clear last axis
for z = j+4:32
    axes(ha(z))
    axis off
end

%save figure
%set(gcf,'Position',[500 1000 900 1150]); %Don't use when plotting all data
%set(gcf,'PaperPositionMode','auto');
set(gcf,'units','normalized','outerposition',[0 0 1 1])
%position = get(gcf,'Position');
%set(gcf,'PaperPosition',[0.5,0,position(3:4)]);
if remove_outliers == 1
    save2pdf([savePath genotypes{i} '_OutliersRemoved.pdf'],gcf)
    %print(gcf,'-dpdf',[folders{1} genotypes{i} '_OutliersRemoved.pdf'])
    saveas(gcf,[savePath genotypes{i} '_OutliersRemoved.fig'])
else
    save2pdf([savePath genotypes{i} '.pdf'],gcf)
    %print(gcf,'-dpdf',[folders{1} genotypes{i} '.pdf'])
    saveas(gcf,[savePath genotypes{i} '.fig'])
end

