function plot_modelsinfolder(folder)
if(~isdeployed)
addpath(genpath('./export_fig'));
addpath(genpath('./tight_subplot'));
end
if strcmp(folder(end),'/') == 0
    folder = [folder '/'];
end

dir_list = dir(folder);
file_num = length(dir_list);

ncolumns = 8;
nrows = 4;
plotnum = 0;
clf;
figure('OuterPosition',[0 0 ncolumns*200 nrows*100]);
ax = tight_subplot(nrows,ncolumns,[.005 .01],[.01 .01],[.01 .01]);

plotnum = 0;

for y = 1:file_num
    file = dir_list(y).name; %pull out the file name
    [~,root,ext] = fileparts(file);
    path_file = [folder file];
    TG = strcmp(ext,'.mat');
    
    
    if TG == 1
            plotnum = plotnum + 1;
            axes(ax(plotnum));
            
            %get plot data and limits
            load(path_file,'pulse_model');
            right = size(pulse_model.fhM,2); 
            bottom_value = min(min(pulse_model.Z2fhM));
            top_value =  max(max(pulse_model.Z2fhM));
            bottom = bottom_value - std(min(pulse_model.Z2fhM));
            top = top_value + std(max(pulse_model.Z2fhM));
            parsed_filename = textscan(root,'%s','Delimiter','_');
            ch = strcat(parsed_filename{1}(2),parsed_filename{1}(3));
            
            if numel(pulse_model.fhM) >0
                %plot data
                hold on
                axis([1 right bottom top])
                plot(pulse_model.Z2fhM','g');
                plot(pulse_model.fhM,'k');
                text(1,bottom+ std(min(pulse_model.Z2fhM)),ch,'FontSize',10);
                hold off
            end
            text(1,bottom+ std(min(pulse_model.Z2fhM)),ch,'FontSize',10);

    end
end
r=regexp(folder,'/','split');
folder_name = char(r(end-1));
outfile = [folder folder_name '_pulsemodels.png'];
warning('off','MATLAB:LargeImage');
export_fig(outfile,'-r300');
