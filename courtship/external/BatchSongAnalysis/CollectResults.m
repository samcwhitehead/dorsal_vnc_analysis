function [ControlResults,Results] = CollectResults(folders,genotypes,control_folders,control_genotypes)

% folders = Cell array of folders containing genotypes to plot
% genotypes = Cell array of genotypes to plot
% control_folders = Cell array of folders containing control
% control_genotypes = Cell array of control genotype(s). All control
% genotypes are included in all plots, but statistics are calculated for
% grand mean of controls.
% remove_outliers = logical (1/0) to plot with outliers removed or original
% data
%
% To work with the output, use [] to grab data from each category
% e.g.
% [Results.SineToNull]
%

%collect control data


numControls = numel(control_genotypes);
numGenotypes = numel(genotypes);

controls = struct();
for i  = 1:numControls
    count = 0;
    for ii = 1:numel(control_folders)
        if control_folders{ii}(end) ~= '/'
            control_folders{ii} = [control_folders{ii} '/'];
        end
        dir_list = dir(control_folders{ii});
        for g = 1:numControls
            idx = ~cellfun('isempty',strfind({dir_list.name},control_genotypes{g}));
            if sum(idx) < 1
                fprintf('Control genotype not found. Check details.\n')
                return
            end
        end
        
        for iii = 1:numel(dir_list)
            file = dir_list(iii).name;
            [~,~,ext] = fileparts(file);
            if ~isempty(strfind(file,control_genotypes{i})) && strcmp(ext,'.mat')
                count = count + 1;
                AR = load([control_folders{ii} file], '-mat');
                if count == 1
                    ControlResults = AR.Analysis_Results;
                else
                    ControlResults = [ControlResults AR.Analysis_Results];
                end
            end
        end
    end
    varname = genvarname(control_genotypes{i});
    controls.(varname) = ControlResults;
end

SampleSize = zeros(1,numControls+1);
for x = 1:numControls
    varname = genvarname(control_genotypes{x});
    SampleSize(x) = numel(controls.(varname));
end


%collect target results

results = struct();
for i  = 1:numel(genotypes)
    count = 0;
    for ii = 1:numel(folders)
        if folders{ii}(end) ~= '/'
            folders{ii} = [folders{ii} '/'];
        end
        dir_list = dir(folders{ii});
        for g = 1:numGenotypes
            idx = ~cellfun('isempty',strfind({dir_list.name},genotypes{g}));
            if sum(idx) < 1
                fprintf('Genotype not found. Check details.\n')
                return
            end
        end
        for iii = 1:numel(dir_list)
            file = dir_list(iii).name;
            [~,~,ext] = fileparts(file);
            if ~isempty(strfind(file,genotypes{i})) && strcmp(ext,'.mat')
                count = count + 1;
                AR = load([folders{ii} file], '-mat');
                if count == 1
                    Results = AR.Analysis_Results;
                else
                    Results = [Results AR.Analysis_Results];
                end
            end
        end
    end
    varname = genvarname(genotypes{i});
    results.(varname) = Results;
end
