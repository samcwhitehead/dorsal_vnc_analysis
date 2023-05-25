function songIndexStruct = getPulseAndSineIndices(data_in, params)
%--------------------------------------------------------------------------
% given an analyzed song, pull out the sine, pulse, and song index that
% tells how often a fly produces each type of song
% (adapted from Arthur code)
%--------------------------------------------------------------------------

% define parameters
minIPI = params.minIPI ; %100;
maxIPI = params.maxIPI ; %3000;
LLR_threshold = params.LLR_threshold ;
pulseCullType = params.pulseCullType ; % 'AmpCull' (more premissive) or 'IPICull' (stricter)
sineCullType = params.sineCullType ; 
pauseThreshold = params.pauseThreshold ; 

% initialize data structure 
songIndexStruct = struct() ; 

%--------------------------------------------------------------------------
%% get pulse data from arthur code analysis
Pulses = data_in.Pulses ;
Data = data_in.Data ; 

switch pulseCullType
    case 'IPICull'
        pulseLLR = [Pulses.Lik_pulse2.LLR_fh] ;
    case 'AmpCull'
        pulseLLR = [Pulses.Lik_pulse.LLR_fh] ;
end

try
    pulses.w0 = Pulses.(pulseCullType).w0(pulseLLR > LLR_threshold );
    pulses.w1 = Pulses.(pulseCullType).w1(pulseLLR > LLR_threshold );
    pulses.wc = pulses.w1 - ...
        ((pulses.w1 - pulses.w0)./2);
    pulses.x = ...
        GetClips(pulses.w0,pulses.w1,Data.d);
catch
    pulses.x = {};
end
%--------------------------------------------------------------------------
%% get sine data from arthur code analysis
try
    sines = data_in.Sines.(sineCullType) ;
    sines.clips = GetClips(sines.start,sines.stop,Data.d)';
catch
    disp('No sine song data!')
    sines.start = [];
    sines.stop = [];
    sines.clips = {};
    keyboard
end

%--------------------------------------------------------------------------
%% calc peak to peak IPIS
try
    p = pulses.wc;
    p_shift_one = circshift(p,[0 -1]);
    ipi.d=p_shift_one(1:end-1)-p(1:end-1);
    ipi.t = p(1:end-1);
    %ipi = fit_ipi_model(pulses);
    %cull IPIs
    culled_ipi.d = ipi.d(ipi.d > minIPI & ipi.d < maxIPI);
    culled_ipi.t = ipi.t(ipi.d > minIPI & ipi.d < maxIPI);
catch
    %ipi.ipi_mean = [];
    %ipi.ipi_SD = [];
    %ipi.ipi_d = [];
    %ipi.ipi_time = [];
    %ipi.fit = {};
    culled_ipi.d = [];
    culled_ipi.t = [];
end

%% calculate song bouts/pauses
if numel(culled_ipi.d) > 1
    %find IPI trains
    IpiTrains = findIpiTrains(culled_ipi);
    %discard IPI trains shorter than max allowed IPI
    IpiTrains.d = IpiTrains.d(cellfun(@(x) ((x(end)-x(1))>maxIPI),IpiTrains.t));
    IpiTrains.t = IpiTrains.t(cellfun(@(x) ((x(end)-x(1))>maxIPI),IpiTrains.t));

    %find All Pauses
    Pauses = findPauses(Data,sines,IpiTrains);
    
    %find Song Bouts
    Bouts = findSongBouts(Data,sines,IpiTrains,Pauses,pauseThreshold);
else
    %IpiTrains.d = {};
    %IpiTrains.t = IpiTrains.d;
    Pauses.PauseDelta = [];
    Pauses.Type = {};
    Pauses.Time = [];
    Pauses.sinesine = [];
    Pauses.sinepulse = [];
    Pauses.pulsesine = [];
    Pauses.pulsepulse = [];
    Bouts.Start = [];
    Bouts.Stop = [];
%     Bouts.x = {};
end

%--------------------------------------------------------------------------
%% Total recording, sine, pulse, bouts
recording_duration = length(Data.d);
if numel(sines.start) > 0
    SineTrainNum = numel(sines.start);
    SineTrainLengths = (sines.stop - sines.start);
    SineTotal = sum(SineTrainLengths);
else
    SineTrainNum = 0;
    SineTrainLengths = 0;
    SineTotal = 0;
end

if numel(IpiTrains.t) > 0
    PulseTrainNum = numel(IpiTrains.t);
    PulseTrainLengths = cellfun(@(x) x(end)-x(1), IpiTrains.t);
    PulseTotal = sum(PulseTrainLengths);
else
    PulseTrainNum = 0;
    PulseTrainLengths = 0;
    PulseTotal = 0;
end

%% calculate song index values 

% across recording values
PulseIndex = PulseTotal / recording_duration ; 
SineIndex = SineTotal / recording_duration ; 
SongIndex = (PulseTotal + SineTotal) / recording_duration ;  
% per unit time values
PulseTrainsPerMin = PulseTrainNum  * 60/(recording_duration / Data.fs);
PulsesPerMin = PulseTotal * 60/(recording_duration / Data.fs);
SineTrainsPerMin = SineTrainNum * 60 /(recording_duration / Data.fs);
SinePerMin = SineTotal * 60/(recording_duration / Data.fs);
BoutsPerMin = numel(Bouts.Start) * 60 / (recording_duration / Data.fs);
SongPerMin = (PulseTotal + SineTotal) * 60/(recording_duration / Data.fs);

%% add data to structure
songIndexStruct.pulseIndex = PulseIndex ; 
songIndexStruct.sineIndex = SineIndex ; 
songIndexStruct.songIndex = SongIndex ; 

songIndexStruct.PulseTrainsPerMin = PulseTrainsPerMin ; 
songIndexStruct.PulsesPerMin = PulsesPerMin ; 
songIndexStruct.SineTrainsPerMin = SineTrainsPerMin ; 
songIndexStruct.SinePerMin = SinePerMin ; 
songIndexStruct.BoutsPerMin = BoutsPerMin ; 
songIndexStruct.SongPerMin = SongPerMin ; 


end