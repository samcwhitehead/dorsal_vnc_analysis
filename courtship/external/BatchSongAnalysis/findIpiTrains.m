function IpiTrains = findIpiTrains(culled_ipi)

%grab indexes for start and stop of bouts (cols 1 and 2 of
%boutIdx)
d = culled_ipi.d;
t = culled_ipi.t;
numIpis = numel(culled_ipi.d);
%set up arrays
startIdx = zeros(numIpis,1);
stopIdx = zeros(numIpis,1);
j = 1;
startIdx(j) = 1;

x = [];
for i = 2:numel(t);
    if t(i) - t(i-1) ~= d(i-1)%collect ipis where the ipi matches the distance between pulses
        stopIdx(j) = i-1;
        j=j+1;
        startIdx(j) = i;
        
        x(j) = t(i);
    end
end
stopIdx(j) = i;

startIdx(startIdx == 0) = [];
stopIdx(stopIdx == 0) = [];

numTrains = numel(startIdx);
Trains = cell(1,numTrains);
Times = cell(1,numTrains);

for i = 1:numTrains
    Trains{i} = d(startIdx(i):stopIdx(i));
    Times{i} = t(startIdx(i):stopIdx(i));
end

IpiTrains.d = Trains;
IpiTrains.t = Times;