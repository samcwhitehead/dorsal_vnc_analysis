function [metadata, summary] = kinefly_analysis(ABF_FILE, desired_trials, protocol_name)

% it would be smarter to find the largest peak or valley in the first few
% seconds after the stimulus instead of looking at a fixed window.
% implement that later.

% CONSTANTS
SAMPLE_RATE=20000; % axoscope acquisition
WBF_THRESHOLD=100; % threshold to use for excluding trials
THRUST_THRESHOLD=0.2; % threshold to use for counting poor fliers
MAX_X_POSITION = 96; % last x position index
DUR_PRESTIM=1*SAMPLE_RATE; %1s
DUR_POSTSTIM=5*SAMPLE_RATE; %5s
XPOS_BINS = linspace(0, 5, MAX_X_POSITION); % divide 0 to 5 volts into bins

% SET UP
close all;
START_PATH = pwd;
[PATHSTR, NAME, EXT] = fileparts(ABF_FILE);

name_split = strsplit(NAME, '_');

metadata.linename = name_split{2};
metadata.experiment_name = NAME;
metadata.experiment_date = name_split{1};

if ~strcmp(PATHSTR,'')
    cd(PATHSTR)
end

[d, ~, ~]=abfload([NAME EXT]);

%% FIND STIMULUS AND STORE ROIS
% finding peaks of channel 16 (stim)
idx_start=find(diff(d(:, 16))>0.05);
% only keep peaks that are 100 indices apart
idx_start=vertcat(idx_start(diff(idx_start)>100), idx_start(end));

% only use the idx_start events for the desired protocol
if nargin > 1
    disp('Only the following trials will be considered:')
    disp(desired_trials)
    idx_start = idx_start(desired_trials);
end

%% EXTRACT RELEVANT TRIAL DATA

for k=1:length(idx_start)
    try
        output{1}(k,:)=d(idx_start(k)-DUR_PRESTIM:idx_start(k)+DUR_POSTSTIM, 1); % LEFT_AMP_C3
        output{2}(k,:)=d(idx_start(k)-DUR_PRESTIM:idx_start(k)+DUR_POSTSTIM, 2); % RIGHT_AMP_C3
        output{3}(k,:)=d(idx_start(k)-DUR_PRESTIM:idx_start(k)+DUR_POSTSTIM, 3); % HEAD_ANG_C3
        output{4}(k,:)=d(idx_start(k)-DUR_PRESTIM:idx_start(k)+DUR_POSTSTIM, 4); % AUX_ROI_C3
        output{5}(k,:)=d(idx_start(k)-DUR_PRESTIM:idx_start(k)+DUR_POSTSTIM, 5); % LEFT_DEV_C2
        output{6}(k,:)=100*d(idx_start(k)-DUR_PRESTIM:idx_start(k)+DUR_POSTSTIM, 6); % LEFT_RAD_C2
        output{7}(k,:)=d(idx_start(k)-DUR_PRESTIM:idx_start(k)+DUR_POSTSTIM, 7); % RIGHT_DEV_C2
        output{8}(k,:)=d(idx_start(k)-DUR_PRESTIM:idx_start(k)+DUR_POSTSTIM, 8); % AUX_ROI_C2
        output{9}(k,:)=d(idx_start(k)-DUR_PRESTIM:idx_start(k)+DUR_POSTSTIM, 9); % RIGHT_DEV_C1
        output{10}(k,:)=100*d(idx_start(k)-DUR_PRESTIM:idx_start(k)+DUR_POSTSTIM, 10); % RIGHT_RAD_C1
        output{11}(k,:)=d(idx_start(k)-DUR_PRESTIM:idx_start(k)+DUR_POSTSTIM, 11); % LEFT_DEV_C1
        output{12}(k,:)=d(idx_start(k)-DUR_PRESTIM:idx_start(k)+DUR_POSTSTIM, 12); % AUX_ROI_C1
        output{13}(k,:)=d(idx_start(k)-DUR_PRESTIM:idx_start(k)+DUR_POSTSTIM, 13); % WBF
        output{14}(k,:)=arrayfun(@(x) find(x < XPOS_BINS, 1, 'first'),...
            d(idx_start(k)-DUR_PRESTIM:idx_start(k)+DUR_POSTSTIM, 14)); % POS_X_IDX
        output{15}(k,:)=d(idx_start(k)-DUR_PRESTIM:idx_start(k)+DUR_POSTSTIM, 15); % POS_Y
        output{16}(k,:)=d(idx_start(k)-DUR_PRESTIM:idx_start(k)+DUR_POSTSTIM, 16); % STIM
        output{17}(k,:)=output{1}(k,:) - output{2}(k,:); % L - R, turn
        output{18}(k,:)=output{1}(k,:) + output{2}(k,:); % L + R, thrust
        output{19}(k,:)=mod(output{14}(k,:), 16); % relative stripe position index
        output{20}(k,:)=pos_hist(output{19}(k,:)); % binned data
        output{21}(k,:)=calc_hwm(output{20}(k,:)); % HWM, histogram width metric for relative stripe position
    catch
    end
end

time_series = ((1:(DUR_PRESTIM+DUR_POSTSTIM+1))/SAMPLE_RATE)';

%% CALCULATE WBF
% % not perfect, but we set a threshold as minimum + 30% range
% min_ = min(min(output{13}(:,15001:25000)));
% max_ = max(max(output{13}(:,15001:25000)));
% threshold_ = min_ + 0.3 * (max_ - min_);
output{13} = find_wbf_spectrogram(d, time_series) ; 

% for k=1:length(output{13}(:,1))
%     output{13}(k,:) = find_wbf_spectrogram(output{13}(k,:), time_series) ; 
% %     idx_event_=[];
% %     idx_event_=find_wbf(output{13}(k,:),threshold_);
% %     if numel(idx_event_) > 2 %1
% %         idx_event_(idx_event_(:,2)<120 | idx_event_(:,2)>280,2) = NaN;
% %         %disp(k)
% %         output{13}(k,:) = interp1(idx_event_(:,1), idx_event_(:,2), time_series); % resampled WBF
% %     else
% %         output{13}(k,:) = zeros(1, length(time_series));
% %     end
% end

%% Throw out events with WBF < threshold
n = 0;
stop_pre = 0;
stop_post = 0;
for k=1:length(idx_start)
    % throw out the trial if the fly wasn't flying in the 1s before or
    % if the stripe wasn't moving, do not increment n
    if isnan(nanmean(nanmean(output{13}(k,1:DUR_PRESTIM)))) || ...
       (nanmean(nanmean(output{13}(k,1:DUR_PRESTIM))) < WBF_THRESHOLD) || ...
       numel(unique(output{14}(k,:))) == 1
        stop_pre = stop_pre + 1;
        continue
    % throw out the trial if the fly spent more than 20% of the poststim
    % duration with WBF=NAN or WBF < WBF_THRESHOLD, or if the fly spent
    % < 20% of the poststim duration with THRUST < THRUST_THRESHOLD
%      elseif (sum(output{13}(k,1:dur_poststim) < WBF_THRESHOLD) + ...
%              sum(isnan(output{13}(k,1:dur_poststim))))> 0.2*dur_poststim || ...
%             sum(output{18}(k,1:dur_poststim) < THRUST_THRESHOLD) > 0.2*dur_poststim
%        stop_post = stop_post + 1;
%        continue
    else
        n = n + 1;
        for i = 1:length(output)
            good_events{i}(n,:) = output{i}(k,:);
        end
    end
end

output = good_events;
metadata.n = n; % number of trials counted
metadata.stop_pre = stop_pre; % number of trials where fly stopped before
metadata.stop_post = stop_post; % number of trials where fly stoped after

%% SAVE CALCULATED MEANS AND STDEV
idx_prestim = (DUR_PRESTIM - (SAMPLE_RATE * 0.5 - 1)):DUR_PRESTIM; % half a second before
idx_stim = (DUR_PRESTIM + 1):(DUR_PRESTIM + 0.1 * SAMPLE_RATE); % duration of stimulus
idx_poststim = (DUR_PRESTIM + 1):(DUR_PRESTIM + (SAMPLE_RATE * 0.5)); % half a second after

for k = 1:19
    summary{k}.pre_val = nanmean(output{k}(:,idx_prestim),2);
    summary{k}.pre_val_mean = nanmean(nanmean(output{k}(:,idx_prestim),2));
    summary{k}.pre_val_std = nanstd(nanmean(output{k}(:,idx_prestim),2));
    summary{k}.stim_val = nanmean(output{k}(:,idx_stim),2);
    summary{k}.stim_val_mean = nanmean(nanmean(output{k}(:,idx_stim),2));
    summary{k}.stim_val_std = nanstd(nanmean(output{k}(:,idx_stim),2));
    summary{k}.post_val = nanmean(output{k}(:,idx_poststim),2);
    summary{k}.post_val_mean = nanmean(nanmean(output{k}(:,idx_poststim),2));
    summary{k}.post_val_std = nanstd(nanmean(output{k}(:,idx_poststim),2));
    summary{k}.stim_delta_mean = nanmean(nanmean(output{k}(:,idx_stim),2)-nanmean(output{k}(:,idx_prestim),2));
    summary{k}.stim_delta_std = nanstd(nanmean(output{k}(:,idx_stim),2)-nanmean(output{k}(:,idx_prestim),2));
    summary{k}.poststim_delta_mean = nanmean(nanmean(output{k}(:,idx_poststim),2)-nanmean(output{k}(:,idx_prestim),2));
    summary{k}.poststim_delta_std = nanstd(nanmean(output{k}(:,idx_poststim),2)-nanmean(output{k}(:,idx_prestim),2));
end

%% SUMMARY PLOTS
if 1
    %%PLOT DIRECTLY MEASURED STATS
    names = {'CAM3 - WBA L', 'CAM3 - WBA R', 'CAM3 - Head', 'CAM3 - Aux.',...
        'CAM2 - L Dev', 'CAM2 - L Rad', 'CAM2 - R Dev', 'CAM2 - Aux.',...
        'CAM1 - R Dev', 'CAM1 - R Rad', 'CAM1 - L Dev', 'CAM1 - Aux.',...
        'WBF', 'X position', 'Y position', 'Stimulus',...
        'CW Torque', 'Thrust', 'dX', 'Stripe Position'};
    units = {'rad', 'rad', 'rad', 'int.',...
        'rad', 'px', 'rad', 'int.',...
        'rad', 'px', 'rad', 'int.',...
        'Hz', 'pos', 'pos', 'V',...
        '', '', 'pix/s', 'pix'};

    figure(1)
    set(1, 'Position', [60, 55, 1200, 800])

    plot_idx = 1;
    for k = [1:2 17:18 3:13]
        % produce mean time series
        mean_val = nanmean(output{k});
        std_val = nanstd(output{k});
        bound_time = time_series(1:200:end)';
        lower = mean_val - std_val;
        lower = lower(1:200:end);
        upper = mean_val + std_val;
        upper = upper(1:200:end);
        lower(isnan(lower)) = 0;
        upper(isnan(upper)) = 0;

        % calculate some useful stats for scaling axes
        min_val = min(mean_val);
        max_val = max(mean_val);
        range_val = range(mean_val);

        %subplot_spaceplots(8,2,plot_idx)
        subplot(8,2,plot_idx)
        hold on

        if k == 20
            for trial = 1:n
                plot(0:15, output{k}(trial,:), 'LineWidth', 0.1)
            end
            ylim([0, 30000])
        else
            % plot the mean time series and set limits
            fill([bound_time,fliplr(bound_time)],[lower,fliplr(upper)],'b','linestyle','none');
            plot(time_series', mean_val, 'b', 'LineWidth', 0.2)
            ylim([(min_val - 0.15 * range_val) (max_val + 0.15 * range_val)])
        end

        xlim([0 time_series(end)]);
        ylabel(units{k})

        yrange = range(ylim);
        ymin = min(ylim);
        ymax = max(ylim);

        % plot means with error bars for 0.5s before and after stimulus
        pre_val = output{k}(:,idx_prestim);
        stim_val = output{k}(:,idx_stim);
        post_val = output{k}(:,idx_poststim);
        errorbar([0.75 1.05 1.25], ...
            [nanmean(nanmean(pre_val,2)) nanmean(nanmean(stim_val,2)) nanmean(nanmean(post_val,2))], ...
            [nanstd(nanmean(pre_val,2)) nanstd(nanmean(stim_val,2)) nanstd(nanmean(post_val,2))], ...
            'ko-', 'LineWidth', 1)

        % plot vertical line at stimulus time
        patch([1 1 1.1 1.1], [ymin ymax ymax ymin], 'r', 'EdgeColor', 'none')
        alpha(0.1)

        % labels and tickmarks
        text(6, ymax, names{k}, 'FontSize', 8, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top')

        % draw pre- and post-stimulus regions on axis
        plot([0.5 1], [ymin ymin], 'r', 'LineWidth', 3)
        plot([1 1.5], [ymin ymin], 'g', 'LineWidth', 3)

        if plot_idx == 14 || plot_idx == 15
            xlabel('time (s)')
            set(gca, 'XTick', linspace(0, 6, 7))
            set(gca, 'YTick', linspace(ymin, ymax, 3))
        else
            set(gca, 'XTick', linspace(0, 6, 7))
            set(gca, 'XTickLabel', [])
            set(gca, 'YTick', linspace(ymin, ymax, 3))
        end

        plot_idx = plot_idx + 1;
    end


    text(8, ymin + yrange*0.5, [NAME ' (N=' num2str(n) ' of ' num2str(length(idx_start)) ' trials)'], ...
        'FontSize', 12, 'Interpreter', 'none')
    text(8, ymin, ['Prestim stop: ' num2str(stop_pre) ', Poststim stop: ' num2str(stop_post)], ...
        'FontSize', 12, 'Interpreter', 'none')

    if nargin > 1
        text(8, ymin - yrange*0.5, protocol_name, ...
            'FontSize', 12, 'Interpreter', 'none')
    end

    % STORE ALL DATA
    saveas(gcf, [NAME '.fig'])
%    tightfig;
    %spaceplots(gcf,[0 0 0 0], [0.02 0.02]);
end

if nargin == 1
    save([NAME '_axo_processed.mat'], 'output', 'summary', 'NAME', 'time_series', 'metadata')
    %save2pdf([NAME '.pdf'], 1, 150)
else
    protocol_name = strrep(protocol_name,' ','_');
    save([NAME '_' protocol_name '_axo_processed.mat'], 'output', 'summary', 'NAME', 'time_series', 'metadata')
    %save2pdf([NAME '_' protocol_name '.pdf'], 1, 150)
end
%close 1

% CLEAN UP
cd(START_PATH)
end

