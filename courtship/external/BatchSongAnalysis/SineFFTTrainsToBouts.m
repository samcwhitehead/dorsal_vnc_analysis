function [time,freq] = SineFFTTrainsToBouts(Bouts,sines,sineMFFT,numTrainsPerBout)

numBouts = numel(Bouts.Start);
CatSineFreq = cell(numBouts,1);
CatSineTime = CatSineFreq;
for i = 1:numBouts
    idx = ismember(sines.start,Bouts.Start(i):Bouts.Stop(i));
    s = find(idx,numTrainsPerBout);
    CatSineFreq{i} = vertcat(sineMFFT.freq{s});
    CatSineTime{i} = horzcat(sineMFFT.time{s});
end

%calculate corr in each sine train
%first, eliminate all trains consisting of fewer than three events
idx = cellfun(@(x) sum(~isnan(x))>2,CatSineTime,'UniformOutput',0);
time = CatSineTime(cell2mat(idx));
freq = CatSineFreq(cell2mat(idx));
