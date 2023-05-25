function [pulsesToAnalyze, pulseMat, pulsesNorm, pulseLabels, h_pclass] = ...
    myPulseClassifier(data_in, params, plotFlag)
%--------------------------------------------------------------------------
% a function to test the murthy lab pulse classifier on VNC lines

% load params
LLR_threshold = params.LLR_threshold ;
max_length = params.max_length ;
total_length = 2* max_length + 1;
pulseCullType = params.pulseCullType ; % 'AmpCull' (more premissive) or 'IPICull' (stricter)

%--------------------------------------------------------------------------
%% get pulse data from arthur code analysis

%Fs = data_in.Data.fs; %Hz
%DataFromStart = data_in.Data.d ;
%T = (1:size(DataFromStart,1))/Fs;
Pulses = data_in.Pulses ;

switch pulseCullType
    case 'IPICull'
        pulseLLR = [Pulses.Lik_pulse2.LLR_fh] ;
    case 'AmpCull'
        pulseLLR = [Pulses.Lik_pulse.LLR_fh] ;
end

pulsesToAnalyze.w0 = Pulses.(pulseCullType).w0(pulseLLR > LLR_threshold );
pulsesToAnalyze.w1 = Pulses.(pulseCullType).w1(pulseLLR > LLR_threshold );
pulsesToAnalyze.wc = pulsesToAnalyze.w1 - ...
    ((pulsesToAnalyze.w1 - pulsesToAnalyze.w0)./2);
pulsesToAnalyze.x = ...
    GetClips(pulsesToAnalyze.w0,pulsesToAnalyze.w1,data_in.Data.d);

%--------------------------------------------------------------------------
%% reshape pulse information into matrix

N_pulses = length(pulsesToAnalyze.x) ;
pulseMat = zeros(N_pulses, total_length) ;

for n = 1:N_pulses
    if ~isempty(pulsesToAnalyze.x{n})
        X = [pulsesToAnalyze.x{n}]';
        T = length(X);
        [~,C] = max(abs(X));%get position of max power
        %flip model is strongest power is negative
        if X(C) <0
            X = -X;
        end
        %center on max power
        left_pad = max_length - C;  %i.e. ((total_length/2) - C)
        right_pad = total_length - T - left_pad;
        if (left_pad < 1) || (right_pad < 1)
            continue
        else
            pulseMat(n,:) = [zeros(1,left_pad) X zeros(1,right_pad)];
        end
    end
end
pulseMat = pulseMat((sum(abs(pulseMat),2)>0),:) ;

%--------------------------------------------------------------------------
%% classify pulses
pulsesNorm = normalizePulses(double(pulseMat)/1000);           % normalize pulses
pulseLabels = classifyPulses(pulsesNorm);       % classify pulses - 0=Pfast, 1=Pslow
fprintf('%d/%d Pslow, %d/%d Pfast.\n', sum(pulseLabels==0), ...
    length(pulseLabels), sum(pulseLabels==1), length(pulseLabels))

%--------------------------------------------------------------------------
%% plot results
if plotFlag
    h_pclass = figure('Name', 'pulse type classification',...
        'PaperPositionMode','auto') ;
    clf
    subplot(311)
    T = (1:size(pulseMat,2))/10;%ms
    plot(T, double(pulseMat')/1000, 'Color', [0 0 0 0.2])
    title('raw pulses')
    subplot(312)
    plot(T, pulsesNorm', 'Color', [0 0 0 0.2])
    title('normalized pulses')
    subplot(313)
    hF = plot(T, pulsesNorm(pulseLabels==1,:)', 'Color', [1 0 0 0.2]);
    hold on
    hS = plot(T, pulsesNorm(pulseLabels==0,:)', 'Color', [0 0 0 0.2]);
    hL = legend([hF(1) hS(1)], {' P_{fast} (red)', 'P_{slow} (black)'},...
        'Box', 'off');
    title([sum(pulseLabels==0) '/' length(pulseLabels) ' Pfast, ' ...
        sum(pulseLabels==1) '/' length(pulseLabels) ' Pslow'])
    axis(gcas, 'tight')
    set(gcas, 'Box', 'off', 'Color', 'none')
else
    h_pclass = nan ; 
end


end
