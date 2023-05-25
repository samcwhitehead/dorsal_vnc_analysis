function sMFFT = findSineMaxFFT(winnowed_sine,Fs)

nfft = 100000;

NumSine = numel(winnowed_sine.clips);
maxFFTFreq = cell(NumSine,1);
maxFFTFreqTime = cell(NumSine,1);
% [poolavail,isOpen] = check_open_pool;
%parfor i = 1:NumSine
    for i = 1:NumSine
    ym = winnowed_sine.clips{i};
    boutStart = winnowed_sine.start(i);
    r = length(ym);
    sec = r/10000;
    if r>1
        if sec < 0.1
            wnd = round(Fs*sec);
            try
                z = resample(ym,Fs,10000);
            catch
                warning('missing')
            end
            [Sn,F] = spectrogram(z,wnd,[],nfft,Fs);
            a = find(F>70 & F<300);
            freq2 = F(a);
            voltage = abs(Sn(a,:));
            
            [~,I] = max(voltage); %I = index of max of the signal between 70-300Hz
            maxFFTFreq{i} = freq2(I); %the frequency with this index
            maxFFTFreqTime{i} = boutStart;
            
        elseif sec > 0.1
            wnd = round(0.1*Fs);
            try
                z = resample(ym,Fs,10000);
            catch
                warning('missing')
            end
            [Sn,F] = spectrogram(z,wnd,[],nfft,Fs);
            a = find(F>70 & F<300);
            freq2 = F(a);
            voltage = abs(Sn(a,:));
            
            [~,I] = max(voltage); %I = index of max of the signal between 70-300Hz
            maxFFTFreq{i} = freq2(I); %the frequency with this index
            maxFFTFreqTime{i} = boutStart:500:boutStart+(500*length(I))-1;
        end
    end
end

sMFFT.freq = maxFFTFreq;
sMFFT.time = maxFFTFreqTime;
sMFFT.freqAll = cell2mat(maxFFTFreq);
sMFFT.timeAll = cell2mat(maxFFTFreqTime');


