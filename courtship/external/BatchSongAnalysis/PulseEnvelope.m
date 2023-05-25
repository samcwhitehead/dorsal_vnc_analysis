function [Start,Center,Stop] = PulseEnvelope(Data, Pulses, MaxExtend, step, NumStd)

%Pulses = most inclusive Pulses structure that includes start and stop
%times
%MaxExtend = maximum extension to start and stop of original pulse to try
%(efault = 150)
%step = step size to try from original pulse to MaxExtend (default  = 10)
%NumStd= # standard deviations of noise to use as cutoff (default = 2)

addpath(genpath('./chronux'))


if nargin < 3
    MaxExtend = 150;
end

if nargin < 4
    step = 10;
end

if nargin < 5
    NumStd  = 2;
end

%fprintf(1,'   Finding Noise Level\n');

if numel(Data.d) > 1e6
    shortdata = Data.d(1:1e6);
else
    shortdata = Data.d;
end
[ssf] = sinesongfinder(shortdata,Data.fs,12,20,.1,.01,.05,[0 1000]); %returns ssf, which is structure containing the following fields
noise = findnoise(ssf,3,80,1000);
StDev = noise.sigma * NumStd;

pulseData = Pulses;

%check here for pulses that start or stop too close to end (within ±
%Extend) -- delete these pulses

KeepIdx = find(pulseData.w0 > MaxExtend & pulseData.w1 < (numel(Data.d) - MaxExtend)); %idx of pulses that are < length of song - MaxExtend
pulseData.w0 = pulseData.w0(KeepIdx);
pulseData.wc = pulseData.wc(KeepIdx);
pulseData.w1 = pulseData.w1(KeepIdx);

AllStarts = pulseData.w0;
NumClips = numel(AllStarts);

for j = 0:step:MaxExtend
    
    Extend = j;
    clips = GetClips(pulseData.w0-Extend,pulseData.w1+Extend,Data.d);
    k = cell2mat(clips);
    SmoothedClips = smooth(k);
    SmoothedClips= reshape(SmoothedClips,size(k,1),size(k,2));
    Env = smooth(abs(hilbert(SmoothedClips)));
    Env = reshape(Env,size(k,1),size(k,2));
    [~,MaxIdx]=max(Env);
    %[~,center] = max(abs(SmoothedClips));
    
    Start = NaN(NumClips,1);
    Stop = Start;
    Center = Start;
    
    for i = 1:NumClips;
        start = find(Env(1:MaxIdx,i) < (StDev),1,'last');
        stop = find(Env(MaxIdx:end,i) < (StDev),1,'first');
        
        %         if strcmp(StartOrCenter,'Start')
        if ~isempty(stop)  &&  ~isempty(start)
            Start(i) = AllStarts(i) - Extend + start; %start at estimated pulse start
            Stop(i) = AllStarts(i) - Extend + stop + MaxIdx(i);
            %             end
            %         elseif strcmp(StartOrCenter,'Center')
            %             if ~isempty(stop)
            Center(i) = pulseData.wc(i); %AllStarts(i) -Extend + center(i); %start at pulse center
            %                 Stop(i) = AllStarts(i) - Extend + stop + MaxIdx(i);
            %             end
        end
    end
    
    Start = Start(~isnan(Start));
    Stop = Stop(~isnan(Stop));
    Center = Center(~isnan(Center));
    if numel(Start) > NumClips * 0.99
        break
    end
end