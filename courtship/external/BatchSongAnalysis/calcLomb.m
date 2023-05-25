function [lombStats] = calcLomb(culled_ipi,fs,alphaThresh)

%
%calc Lomb periodgram sign
%peaks at certain alpha (use 0.01 to start)
%


%calculate lomb-scargle periodgram
[P,f,alpha]=lomb(culled_ipi.d,culled_ipi.t./fs);
%get peaks
peaks = regionalmax(P);
%get f,alpha,Peaks for peaks < desired alpha
fPeaks = f(peaks);
alphaPeaks = alpha(peaks);

signF = fPeaks(alphaPeaks < alphaThresh);
signAlpha = alphaPeaks(alphaPeaks <alphaThresh);
signPeaks = P(alphaPeaks < alphaThresh);

lombStats.F = signF;
lombStats.Alpha = signAlpha;
lombStats.Peaks = signPeaks;
