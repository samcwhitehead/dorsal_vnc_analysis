function Plot_PulseSegmentation_Results(data,winnowed_sine,pulseInfo,pulseInfo2)
%USAGE Plot_PulseSegmentation_Results(data,winnowed_sine,pulseInfo,pulseInfo2)
%This is the primary utility for examining the direct output of ProcessSong. Does not plot likelihoods.
figure; clf;
%plot the signal
plot((1:size(data.d,1))./data.fs,data.d,'Color',[.742 .742 .742])
hold on; 

%plot sine data
if ~isempty(winnowed_sine) && numel(winnowed_sine.start) > 0 
    if winnowed_sine.start>0
        for n = 1:size(winnowed_sine.start,1)
            x_start = round(winnowed_sine.start(n));
            x_stop = round(x_start + size(winnowed_sine.clips{n},1));
            time = (x_start:x_stop-1);
            y = winnowed_sine.clips{n};
            plot(time./data.fs,y,'b')
        end
    end
end
zoom xon;


if numel(pulseInfo) > 0
    xpls = nan*data.d; %vector of NaNs the length of xs
    for i = 1:length(pulseInfo2.x);
        a = pulseInfo2.w0(i);
        b = pulseInfo2.w1(i);
        t = (a:b);
        y = data.d(a:b);
        plot(t./data.fs,y,'r'); %hold on;
    end
    
%plot pulse data
n = length(pulseInfo2.x);
wc = nan(1,n); %NaNs the length of n (number of pulses in pulseInfo2)
mxv = nan(1,n);
mmxv = nan(1,n);
for i = 1:n
        wc(i) = pulseInfo2.wc(i)/data.fs;
        wwc(i) = pulseInfo.wc(i)/data.fs; %pulses before second round of winnowing
        mxv(i) = max(abs(pulseInfo2.x{i}));
        mmxv(i) = max(abs(pulseInfo.x{i})); %pulses before second round of winnowing
end
plot(wwc,mmxv,'.m');
hold on;
plot(wc,mxv,'k^','MarkerFaceColor','k');
end

title('Sine & Pulses','FontSize',14);
hold off
zoom xon;
title('Putative Pulse Regions','FontSize',14);

