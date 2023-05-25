%--------------------------------------------------------------------------
% script to retrieve and store data for plotting summary stats of kinefly
% data
%--------------------------------------------------------------------------
[mfilePath, ~, ~] = fileparts(mfilename('fullpath')) ; 
parentDirectory = fileparts(mfilePath) ;
dataPath = fullfile(parentDirectory, 'data') ;
savePath = dataPath ;

saveFlag = true ;
excludeFlag = true ;
diff_cut_off = 4.0 ; % threshold for derviative in absolute units
diff_cut_off_z = 3.0 ; % 3.0 % threshold for derviative in std units

kf_data_struct = struct() ;
define_constants_kf ;
%--------------------------------------------------------------------------
%% timing info
global_t = T_START:(1/SAMPLE_RATE):T_END ;
sub_time = global_t(1:N_SUB_SAMPLE:end);
t_1_ind = find(sub_time == LED_start_time) ;

plot_idx = (sub_time >= 0.75) & (sub_time <= 1.5) ;
stim_idx = (sub_time >= 1.0) & (sub_time <= 1.25) ; 
%--------------------------------------------------------------------------
%% loop through data directory to get results for each genotype
dataDir = dir([rootPath 'SS*']) ;
driver_cell = arrayfun(@(x) x.name, dataDir,'UniformOutput',0) ;

for i = 1:length(driver_cell)
    %-------------------------------
    % load data for current driver
    driver = driver_cell{i} ;
    dataFolder = fullfile(dataDir(i).folder, dataDir(i).name) ;
    grand_mean_struct = importdata(fullfile(dataFolder, ...
        'grand_mean_struct.mat')) ;
    
    kf_data_struct(i).driver = driver ; 
%     if strcmp(driver,'SS48311') 
%         keyboard ;
%     end
    %----------------------------------------------
    %% exclude bad trials?
    bad_ind_cell = cell(length(grand_mean_struct), 1) ;
    if excludeFlag
        % loop through variables in grand_mean_struct
        for j = 1:length(grand_mean_struct)
            %var_curr = grand_mean_struct(j).var_name ;
            mean_mat = grand_mean_struct(j).data_means ;
            grand_mean = grand_mean_struct(j).grand_mean ;
            
%             diff_mat = bsxfun(@minus, mean_mat(:, plot_idx), ...
%                 grand_mean(plot_idx)) ;
            diff_mat = diff(mean_mat(:,plot_idx),1,2) ;
            diff_mat_norm = nanmax(abs(diff_mat),[],2) ; % L1 norm
            diff_mat_zScore = (diff_mat_norm - nanmean(diff_mat_norm)) ./ ...
                nanstd(diff_mat_norm) ;
            
            bad_ind = (abs(diff_mat_zScore) >= diff_cut_off_z) | ...
                (abs(diff_mat_norm) >= diff_cut_off) ; 
            if (0)
                figure ;
                plot(diff_mat_zScore, 'ko-')
                keyboard
            end
            % hack fix for now, need to update for real
            if (j ==7) && (i == 13)
               bad_ind(7) = true ;  
            end
            bad_ind_cell{j} = bad_ind ; 
        end
    else
        bad_ind_cell{j} = true(size(grand_mean_struct(1).data_means, 1), 1) ; 
    end
    
    % aggregate bad indices
    bad_ind_curr = false(size(grand_mean_struct(1).data_means, 1), 1) ; 
    for p = 1:length(bad_ind_cell)
       bad_ind_curr = (bad_ind_curr | bad_ind_cell{p}) ;  
    end
    if sum(~bad_ind_curr) < 3
        keyboard
    end
    % add to data structure
    kf_data_struct(i).bad_ind = bad_ind_curr ; 
    
    %----------------------------------------------------------------------
    %% loop through data to get scalar stats
    for k = 1:length(grand_mean_struct)
        var_curr = grand_mean_struct(k).var_name ;
        mean_mat = grand_mean_struct(k).data_means ;
        %prestim_mean = grand_mean_struct(k).pre_stim_mean ;
        
        % store data matrix in structure 
        kf_data_struct(i).([var_curr '_mean_mat']) = mean_mat ; 
        %kf_data_struct(i).([var_curr '_pre_stim_mean']) = prestim_mean ; 
        
        % exclude bad trials
        mean_mat = mean_mat(~bad_ind_curr,:) ; 
        
        % get peak change for each trial
        peak_mat = nan(size(mean_mat,1),1) ; 
        for m = 1:size(mean_mat,1)
           data_curr = mean_mat(m,:) ; 
           [max_val, max_ind] = max(abs(data_curr(stim_idx))) ;
           max_val_sign = sign(data_curr(max_ind + t_1_ind)) ; 
           max_val = max_val * max_val_sign ; 
           
           peak_mat(m) = max_val ; 
        end
        
        % get mean and standard deviation of peaks
        mean_peak = nanmean(peak_mat) ; 
        std_peak = nanstd(peak_mat) ; 
        
        % add data to structure
        kf_data_struct(i).([var_curr '_peak_mat']) = peak_mat ; 
        kf_data_struct(i).([var_curr '_peak_mean']) = mean_peak ; 
        kf_data_struct(i).([var_curr '_peak_std']) = std_peak ; 
    end
   fprintf('Completed analysis of %s \n', driver)
    
end
%--------------------------------------------------------------------------
%% save results?
if saveFlag
   save(fullfile(savePath,'kf_data_struct.mat'),'kf_data_struct') 
end