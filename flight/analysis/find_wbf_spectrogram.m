%--------------------------------------------------------------------------
% function to calculate wingbeat frequency based on spectrogram analysis
% rather than the method used in the standard kinefly_analysis.m
%
% INPUTS:
%   - x: the estimated wbf signal from kinefly_analysis.m (output{13})
%   - times_series: time values corresponding to measurements
%   - debugFlag: boolean, should we show plots or not?
%   
%--------------------------------------------------------------------------
function wbf_out = find_wbf_spectrogram(wbf_in, time_series, debugFlag)
%--------------------------------------------------------------------------
%% INPUTS
if ~exist('debugFlag','var')
    debugFlag = false ; 
end
%-------------------------------------------------------------------------- 
%% CONSTANTS
%----------------------
% general params
SAMPLE_RATE=20000; % axoscope acquisition
WBF_THRESHOLD=120; % threshold to use for excluding trials

%-------------------------
% params for spectrogram
freqRange = [100, 300] ; % Hz
freqRes = 30 ; % frequency resolution, in Hz (ends up being finer than this)
overlapPerc = 90 ; % amount of spectrogram overlap

%----------------------
% params for filter
fs = SAMPLE_RATE ; 
f0 = 60 ; % first fundamental of frequency we want to filter out
filterOrder = floor(fs/f0) ; % number of comb notches
Q = 35 ; % quality factor. important? who knows...

BW = (f0/(fs/2))/Q; % filter bandwidth
[filt_b, filt_a] = iircomb(filterOrder, BW,'notch') ; 

% % -------------------------------------------------------------------------
% %% FIRST DO THE OPERATION THAT THE ORIGINAL KINEFLY ANALYSIS DID
% % NB: this doesn't really affect the results, but we want exact consistency
% % set a threshold as minimum + 30% range
% min_ = min(min(x(:,15001:25000)));
% max_ = max(max(x(:,15001:25000)));
% threshold_ = min_ + 0.3 * (max_ - min_);
% 
% idx_event_=find_wbf(x(k,:),threshold_);
% if numel(idx_event_) > 2 %1
%     idx_event_(idx_event_(:,2)<120 | idx_event_(:,2)>280,2) = NaN;
%     %disp(k)
%     x(k,:) = interp1(idx_event_(:,1), idx_event_(:,2), time_series); % resampled WBF
% else
%     x(k,:) = zeros(1, length(time_series));
% end

%--------------------------------------------------------------------------
%% FILTER WBF ESTIMATE AND CALCULATE SPECTROGRAM
try
    % use notch filter to remove 60 Hz noise from current 
    wbf_filt = filtfilt(filt_b, filt_a, wbf_in) ;
    % get cpectrogram
    [p, f, t] = pspectrum(wbf_filt,SAMPLE_RATE,'spectrogram',...
        'FrequencyLimits',freqRange, 'FrequencyResolution',freqRes,...
        'OverlapPercent',overlapPerc) ;
    % find maximum power denisty for each time bin
    [~, f_max_ind] = nanmax(p) ;
    
    % plot results?
    if debugFlag
        figure ;
        imagesc(t,f,p)
        set(gca,'ydir','normal')
        hold on
        plot(t, f(f_max_ind), 'r-')
    end
    
    % get estimate for frequency and interpolate
    wbf_est = f(f_max_ind) ;
    wbf_left_idx = find(~isnan(wbf_est),1,'first') ;
    wbf_left = wbf_est(wbf_left_idx) ;
    wbf_right_idx = find(~isnan(wbf_est),1,'last') ;
    wbf_right = wbf_est(wbf_right_idx) ;
    
    interp_t = [time_series(1) ; t ; time_series(end)] ;
    wbf_pad = [wbf_left ; wbf_est ; wbf_right] ;
    wbf_interp = interp1(interp_t, wbf_pad, time_series, 'spline') ;
    
    % remove bad points
    bad_idx = (wbf_interp < WBF_THRESHOLD) ;
    wbf_interp(bad_idx) = nan ;
    
    wbf_out = wbf_interp ;
    
catch
    wbf_out = [] ;
end

end