function TranProbMatrix = TranProb(Data,sines,Pulses)
fs = Data.fs;
minPause = 0.5*fs;

%collect known sine and pulse events
%pauses == 1, sines == 2, pulses == 3
X = ones(1,numel(Data.d));
for i = 1 : numel(sines.start)
    X(sines.start(i):sines.stop(i)) = 2;
end
for ii = 1:numel(Pulses.w0)
    X(Pulses.w0(ii):Pulses.w1(ii)) = 3;
end

%convert pauses < minPause to previous neighbor (maybe later write to split
%difference)
x = SplitVec(X);
for iii = 2:numel(x)-1
    if x{iii}(1) == 1 %if a pause
        if numel(x{iii}) < minPause %if pause shorter than minPause
            previousState = x{iii-1}(1);
            %followingState = x{iii+1}(1);
            y = zeros(1,numel(x{iii}));
            y(:) = previousState;
            cat(2,x{iii-1},y);
            x{iii} = [];
            
        end
    end
end
x(cellfun('isempty',x)) = [];
%z = cellfun(@(x) downsample(x,1e2),x,'UniformOutput',0);
z = cell2mat(x);
zz = downsample(z,1e2);

%sparse trick to get tran prob (matlabcentral # 41324
t = sparse(zz(1:end-1),zz(2:end),1);
tt = full(t);
e = tt ./repmat(sum(tt,2),1,size(tt,2));

TranProbMatrix = e;

% TRGUESS = [ 0.9162    0.0138    0.0700
%     0.0000    1.0000    0.0000
%     0.0000    0.0000    1.0000];
% EMITGUESS = [   1 0 0
%     0 1 0
%     0 0 1];
% [estTransProb,estEmiss] = hmmtrain(zz,TRGUESS,EMITGUESS);
