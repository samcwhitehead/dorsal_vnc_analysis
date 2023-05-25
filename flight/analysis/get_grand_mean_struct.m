%--------------------------------------------------------------------------
% function that outputs a structure containing the (subsampled) data for a
% given set of data folders, along with the per-fly mean and the grand mean
% across data
%--------------------------------------------------------------------------
function grand_mean_struct = get_grand_mean_struct(dataPaths)
%--------------------------------------------------------------------------
%% params
define_constants_kf

var_str = {'R_AMP', 'L_AMP', 'R_DEV_1', 'L_DEV_2',...
    'R_DEV_2', 'L_DEV_1', 'WBF'} ;
angle_conversion_flag = [1, 1, 1, 1, 1, 1, 0] ;
N_var = length(var_str) ;

% time info
global_t = T_START:(1/SAMPLE_RATE):T_END ;
sub_time = global_t(1:N_SUB_SAMPLE:end);
t_1_ind = find(sub_time == LED_start_time) ;

%initialize struct
grand_mean_struct = struct() ;
cc = 1 ; % counter for each experiment
%% loop through data paths
%--------------------------------------------------------------------------
for i = 1:length(dataPaths)
    dataPath = dataPaths{i} ;
    dataPathSplit = strsplit(dataPath,'\') ;
    %condition = dataPathSplit{end} ;
    
    dirTemp = dir([dataPath '\time_series_struct*']) ;
    if length(dirTemp) == 1
        time_series_struct = importdata(fullfile(dirTemp(1).folder, ...
            dirTemp(1).name)) ;
    else
        disp('No time series structure to load')
        continue
    end
    
    % check if time_series_struct is empty
    if ~isfield(time_series_struct, var_str{1})
        disp('empty time series struct')
        continue
    end
    for k = 1:N_var
        data_mean_mat = nan(length(time_series_struct),length(sub_time)) ;
        data_mat_all = [] ;
        data_ind = [] ;
        pre_stim_mean_mat = [] ; 
        for p = 1:length(time_series_struct)
            if angle_conversion_flag(k)
                data_temp = RAD2DEG*time_series_struct(p).(var_str{k}) ;
            else
                data_temp = time_series_struct(p).(var_str{k}) ;
            end
            
            % if the data we're looking at is WBF, need to do some extra
            % processing for extrapolation errors, etc.
            if (k == N_var)
                data_temp_diff = diff(data_temp,1,2) ;
                
                bad_ind_diff = [true(size(data_temp,1),1), ...
                    (abs(data_temp_diff) >= FREQ_DIFF_THRESH)] ;
                bad_ind_range = (data_temp < FREQ_RANGE(1)) | ...
                    (data_temp > FREQ_RANGE(2)) ;
                bad_ind = bad_ind_diff | bad_ind_range ;
                data_temp(bad_ind) = nan ;
            end
            % apply median filter
            if (k ~= N_var)
                time_series_mat_mf = medfilt1(data_temp,MEDFILT_WINDOW,[],2,...
                    'omitnan','truncate') ;
            else
                time_series_mat_mf = medfilt1(data_temp,301,[],2,...
                    'omitnan','truncate') ;
            end
            % sub sample (alternative)
            time_series_mat_mf = time_series_mat_mf(:,1:N_SUB_SAMPLE:end) ;
            
            % subtract mean from each time series
            pre_stim_mean = nanmean(time_series_mat_mf(:,1:t_1_ind),2) ;
            
            %pre_stim_std = nanstd(time_series_mat_mf(:,1:t_1_ind),2) ;
            time_series_mat_mf = time_series_mat_mf - repmat(pre_stim_mean,1,...
                size(time_series_mat_mf,2)) ;
            
            % interpolate over nans
            if sum(isnan(time_series_mat_mf(:))) > 0
                nan_rows = find(any(isnan(time_series_mat_mf),2)) ;
                for nr = 1:length(nan_rows)
                    nan_ind = isnan(time_series_mat_mf(nr,:)) ;
                    try
                        time_series_mat_mf(nr,nan_ind) = interp1(sub_time(~nan_ind),...
                            time_series_mat_mf(nr,~nan_ind),sub_time(nan_ind),...
                            'linear',0) ;
                    catch
                        time_series_mat_mf(nr,:) = nan ;
                    end
                end
            end
            
            data_mat_all = [data_mat_all ; time_series_mat_mf] ;
            data_mean_mat(p,:) = nanmean(time_series_mat_mf,1) ;
            data_ind = [data_ind ; cc*ones(size(time_series_mat_mf,1),1)] ;
            pre_stim_mean_mat = [pre_stim_mean_mat ; pre_stim_mean] ; 
            cc = cc + 1 ;
        end
        
        if i == 1
            grand_mean_struct(k).var_name = var_str{k} ;
            grand_mean_struct(k).data_mat = data_mat_all ;
            grand_mean_struct(k).data_ind = data_ind ;
            grand_mean_struct(k).data_means = data_mean_mat ;
            grand_mean_struct(k).pre_stim_mean = pre_stim_mean_mat ;
        else
            grand_mean_struct(k).data_mat = [grand_mean_struct(k).data_mat ; ...
                data_mat_all] ;
            grand_mean_struct(k).data_ind = [grand_mean_struct(k).data_ind ; ...
                data_ind ];
            grand_mean_struct(k).data_means = ...
                [grand_mean_struct(k).data_means ; data_mean_mat] ;
            grand_mean_struct(k).pre_stim_mean = ...
                [grand_mean_struct(k).pre_stim_mean ; pre_stim_mean_mat] ;
        end
        
        % test to make sure the variables aren't insane
        if (0)
            figure ;
            hold on
            for q = 1:size(data_mean_mat,1)
                plot(sub_time, grand_mean_struct(k).data_means(q,:))
            end
        end
        
    end
    fprintf('Completed %d / %d data folders \n',i,length(dataPaths))
end

%% calculate grand mean and confidence intervals
% bad_ind_curr = bad_mean_ind_cell{rp} ;
% if ~isempty(bad_ind_curr)
%     good_ind = true(size(data_mean_mat,1),1) ;
%     good_ind(bad_ind_curr) = false ;
%     CI_grand_mean = bootci(N_BOOT_SAMPLES,@nanmean,data_mean_mat(good_ind,:)) ;
%     grand_mean = nanmean(data_mean_mat(good_ind,:)) ;
% else
%     CI_grand_mean = bootci(N_BOOT_SAMPLES,@nanmean,data_mean_mat) ;
%     grand_mean = nanmean(data_mean_mat) ;
% end
disp('Calculating mean and bootstrap confidence intervals...')
for kk = 1:N_var
    % get bootstrapped mean and ci across all samples
    CI_mean_all = bootci(N_BOOT_SAMPLES,@nanmean, ...
        grand_mean_struct(kk).data_mat) ;
    bootstrp_mean = bootstrp(N_BOOT_SAMPLES,@nanmean, ...
        grand_mean_struct(kk).data_mat) ;
    
    % get grand mean and bootstrapped CI
    CI_grand_mean = bootci(N_BOOT_SAMPLES,@nanmean, ...
        grand_mean_struct(kk).data_means) ;
    grand_mean = nanmean(grand_mean_struct(kk).data_means) ;
    % make sure to store this in struct
    grand_mean_struct(kk).grand_mean = grand_mean ;
    grand_mean_struct(kk).mean_all_bootstrp = nanmean(bootstrp_mean) ;
    grand_mean_struct(kk).CI_grand_mean = CI_grand_mean ;
    grand_mean_struct(kk).CI_mean_all = CI_mean_all ;
end
disp('Done')

end