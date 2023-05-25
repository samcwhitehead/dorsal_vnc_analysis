% -------------------------------------------------------------------------
% Script to execute 'kinefly_analysis.m' on all data folders within
% dataroot. This process takes .abf files as input and outputs .mat files
% containing the processed result
%
% -------------------------------------------------------------------------
% Root directory containing all experimental data
[mfilePath, ~, ~] = fileparts(mfilename('fullpath')) ; 
parentDirectory = fileparts(mfilePath) ;
dataroot = fullfile(parentDirectory, 'data') ;
START_PATH = dataroot;

% Folder containing data for a given experiment type
genotype = 'SS01062' ;  % 'SS01062' | 'SS37246' | ...
exptype = 'closed_loop_10x';
exproot = fullfile(dataroot, genotype, exptype, filesep);

% Loop through each experiment
if ~strcmp(exptype,'open_loop')
    cd(exproot);
    folders = dir('20*');
    k = 1
    for idx = 1:length(folders)
        folder = fullfile(exproot, folders(idx).name);
        cd(folder);
        abf_file = dir('*.abf');
        try
            disp(abf_file.name)
            abf_filename_split = strsplit(abf_file.name,'.') ; 
            abf_file_prefix = abf_filename_split{1} ; 
            if (exist([abf_file_prefix '_axo_processed.mat'],'file') == 2)
               disp('Already analyzed; skipping')
               continue ; 
            end
            % Analyze experiment, and store metadata and summary level
            % information
            [data{k}.metadata, data{k}.summary] = kinefly_analysis(abf_file.name);
            k = k+1
        catch 
            warning('unable to process video');
        end
    end
    cd(START_PATH);
    %save([datestr(now,30) '_' exptype '.mat'], 'data');
    clear data
else
    for prot_idx = 1:4
        cd(exproot);
        folders = dir('20*');
        k = 1
        for idx = 1:length(folders)
            folder = fullfile(exproot, folders(idx).name);
            cd(folder);
            name_split = strsplit(folders(idx).name, '_');
            linename = name_split{2};
            abf_file = dir('*.abf');
            mat_file = dir(['*' linename '.mat']);
            try
                disp(abf_file.name)
                disp(mat_file.name)
                load(mat_file.name)
                trial_idx = protocol_order == prot_idx;
                protocol_name = protocol(prot_idx).name;
                [data{k}.metadata, data{k}.summary] = kinefly_analysis(abf_file.name, trial_idx, protocol_name)
                k = k+1
            catch
                disp('failed!')
            end
        end
        cd(START_PATH);
        save([datestr(now,'yyyymmddTHHMMSS') '_' exptype '_' protocol_name '.mat'], 'data')
        clear data
    end
end
