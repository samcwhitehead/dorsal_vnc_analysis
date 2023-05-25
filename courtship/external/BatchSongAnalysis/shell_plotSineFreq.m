if isempty(Sines.LengthCull.clips)
    Sines.LengthCull.clips = GetClips(Sines.LengthCull.start,Sines.LengthCull.stop,Data.d);
end
if isempty(Pulses.AmpCull.x)
    Pulses.AmpCull.x = GetClips(Pulses.AmpCull.w0,Pulses.AmpCull.w1,Data.d);
end

figure(1)
plot((1:length(Data.d))./Data.fs,Data.d,'Color',[.742 .742 .742]);
hold on;
num_events=length(Sines.LengthCull.start);
if num_events>0
    for n = 1:size(Sines.LengthCull.start,1)
        x_start = round(Sines.LengthCull.start(n));
        x_stop = round(x_start + size(Sines.LengthCull.clips{n},1));
        time = (x_start:x_stop-1);
        y = Sines.LengthCull.clips{n};
        plot(time./Data.fs,y,'b')
    end
end



if numel(Pulses.AmpCull) > 0
    for i = 1:length(Pulses.AmpCull.x);
        a = Pulses.AmpCull.w0(i);
        b = Pulses.AmpCull.w1(i);
        t = (a:b);
        y = Data.d(a:b);
        plot(t./Data.fs,y,'r'); %hold on;
    end

end

plot(sineMFFT.timeAll*1e4,sineMFFT.freqAll./100,'.k')

