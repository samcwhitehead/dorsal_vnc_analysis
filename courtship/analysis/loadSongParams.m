%--------------------------------------------------------------------------
% function to load parameters for song analysis and (ideally) keep things
% consistent across different types of analysis 
%
% To do:
%   -right now i just have one set of parameters. but if i want to change
%   this, i should define multiple sets, and then make sure that any time I
%   save something I'm not overwriting analysis with different parameters
%--------------------------------------------------------------------------
function params = loadSongParams(paramCase) 

if nargin < 1
    paramCase = 1 ;
end

switch paramCase
    case 1
        %% which types of sine and pulse song to take
        
        pulseCullType = 'IPICull' ; % either 'AmpCull' or 'IPICull'
        sineCullType = 'LengthCull' ; % either 'LengthCull' or 'PulseCull'
        LLR_threshold = 50 ; %log likelihood ratio used for pulse identification
        
        %% pulse classifier 
        max_length = 125 ; % size to use for normalized pulse matrix
        
        %% song bout detection
        minIPI = 100 ; % minimum inter-pulse interval for song bout (in idx)
        maxIPI = 3000 ; % maximum
        pauseThreshold = 0.5e4; %minimum pause between bouts
        
        %% sine frequency analysis
        f_min = 70 ; % minimum sine song frequency (Hz)
        f_max = 300 ; % maximum
        windowSize = 0.1 ; % spectrogram hamming window size (seconds)
        max_sine_length = 1200 ; % longest sine train to consider for averaged analysis (ms)
        sine_t_bin_size = 1000*windowSize/2 ; % bin size used for spectrogram (ms)
        maxBoutNum = 8 ; % max number of consecutive sine trains to analyze per bout
        minSineInBout = 3 ; % min number of sine trains in bout to consider for analysis
        
        %% plotting preferences (needs to be filled out and fully implemented)
        expr_color = [177,0,38]/255 ;
        ctrl_color = [3,78,123]/255 ;
        
    otherwise
        keyboard
end

%% load parameters to structure
params = struct() ;

params.pulseCullType = pulseCullType;
params.sineCullType = sineCullType ; 
params.LLR_threshold = LLR_threshold ;

params.max_length = max_length ; 

params.minIPI = minIPI ; 
params.maxIPI = maxIPI ;
params.pauseThreshold = pauseThreshold ; 

params.f_min = f_min ;
params.f_max = f_max ;
params.windowSize = windowSize ; 
params.max_sine_length = max_sine_length ;
params.sine_t_bin_size = sine_t_bin_size ;
params.maxBoutNum = maxBoutNum ;
params.minSineInBout = minSineInBout ;

params.expr_color = expr_color ; 
params.ctrl_color = ctrl_color ; 

end