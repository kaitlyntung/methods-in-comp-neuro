%% Load data and set variables
load('Data\retinaData.mat');

%% 3.1 Spike Triggered Average
% Part A

% data cleaning (remove 0 entries)
neuron = 1;
idx = retinaData.stimulusFrameTimes ~= 0;
stimulusFrameTimes = retinaData.stimulusFrameTimes(idx);
stimulus = retinaData.stimulusFrames(:,:,idx);
stimulusTimes = stimulusFrameTimes * 0.1; 
spikeTimes = retinaData.spikes{neuron};

lags = -10:10:300;
strf = zeros(40, 40, length(lags));
spikeCount = 0;

for i = 1:length(spikeTimes)
    spike = spikeTimes(i);
    for lag = 1:length(lags)
        time = spike - lags(lag);
        [~, idx] = min(abs(stimulusTimes - double(time)));
        if ~isempty(idx)
            strf(:, :, lag) = strf(:, :, lag) + stimulus(:, :, idx);
        end
    end
    spikeCount = spikeCount + 1;
end
strf = strf / spikeCount;

figure;
max_val = max(abs(strf(:)));
colormap(nebula(256));
for k = 1:32
    img = strf(:,:,k);
    img_scaled = (img + max_val) / (2 * max_val);
    image(img_scaled);
    axis off;
end