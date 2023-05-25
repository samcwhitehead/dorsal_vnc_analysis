function pMFFT = findPulseModelMaxFFT(pulseModel,Fs)

minFreq = 100;
maxFreq = 1500;
%rewrite in form
L = numel(pulseModel);
NFFT = Fs;
pFFT = fft(pulseModel,NFFT)/L;
%reduce each entry to NFFT/2+1
V = 2*abs(pFFT(1:NFFT/2+1));
[~,i] = max(smooth(V(minFreq:maxFreq)));%exclude low and high freqs
i = i + minFreq;%add back low freq range
pMFFT =i;
