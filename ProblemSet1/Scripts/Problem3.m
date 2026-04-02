%% Load data and set variables
load('Data\retinaData.mat');

%% 3.1 Spike Triggered Average
% Part A
function strf = plot_strf(retinaData, neuron)
% Data cleaning
valid = retinaData.stimulusFrameTimes ~= 0;
stimulusFrameTimes = retinaData.stimulusFrameTimes(valid);
stimulus = double(retinaData.stimulusFrames(:,:,valid)) - 0.5;
stimulusTimes = double(stimulusFrameTimes) * 0.1;
spikeTimes = double(retinaData.spikes{neuron}) * 0.1;

lags = 300:-10:-10;
nLags = length(lags);
spikeCount = length(spikeTimes);

% Pre-flatten stimulus: (1600 x nFrames)
stimFlat = reshape(stimulus, 1600, []);

% Find nearest frame index for every spike ONCE
spikeFrameIdx = round(interp1(stimulusTimes, 1:length(stimulusTimes), ...
                               spikeTimes, 'nearest', 'extrap'));
spikeFrameIdx = max(1, min(spikeFrameIdx, length(stimulusTimes)));

strf = zeros(40, 40, nLags);
for lag = 1:nLags
    lagFrames = round(lags(lag) / (1000/30));
    offsetIdx = spikeFrameIdx - lagFrames;
    valid_idx = offsetIdx >= 1 & offsetIdx <= size(stimFlat, 2);
    strf(:,:,lag) = reshape(sum(stimFlat(:, offsetIdx(valid_idx)), 2), 40, 40);
end
strf = strf / spikeCount;

% Plotting
figure;
max_val = max(abs(strf(:)));
for k = 1:nLags
    subplot(6, 6, k);
    img_scaled = (strf(:,:,k) + max_val) / (2 * max_val);
    image(img_scaled * 255);
    axis image off;
    title(sprintf('%d ms', lags(k)), 'FontSize', 7);
end
colormap(gray(256));
sgtitle(sprintf('STRF - Neuron %d', neuron));

end

%Part B
neurons = [1, 4, 11, 15, 26, 51, 80, 84, 96, 105];
strfs = cell(length(neurons), 1);

for i = 1:length(neurons)
    strfs{i} = plot_strf(retinaData, neurons(i));
end

%% 3.2 Receptive Field Models
% Part A
neurons = [1, 4, 11, 15, 26, 51, 80, 84, 96, 105];
targetLag = 80; % ms — index to plot in Part B
lagIdx80 = find(lags == targetLag);

% Storage for fitted parameters
allParams = zeros(length(neurons), 4 + 32); % [x0, y0, sigma, offset, A_1...A_32]
allFitted = cell(length(neurons), 1);

for n = 1:length(neurons)
    neuron = neurons(n);
    strf = plot_strf(retinaData, neuron);  % your existing function
    
    % --- Build xdata: pixel coordinate grid, replicated for each lag ---
    [X, Y] = meshgrid(1:40, 1:40);
    xdata = [X(:), Y(:)];  % (1600 x 2) grid, same for all lags
    
    % --- ydata: flatten STRF to vector (1600*32 x 1) ---
    ydata = strf(:);  % MATLAB flattens column-major: pixel, then lag
    
    % --- Model function: time-varying 2D Gaussian ---
    % params = [x0, y0, sigma, A_1, A_2, ..., A_32]
    % shared: x0, y0, sigma | free: one amplitude per lag
    gaussModel = @(params, xdata) reshape( ...
        cell2mat(arrayfun(@(k) ...
            params(3+k) * exp(-((xdata(:,1)-params(1)).^2 + ...
                               (xdata(:,2)-params(2)).^2) / (2*params(3)^2)), ...
            1:32, 'UniformOutput', false)'), [], 1);
    
    % --- Warm start: initialize from the lag with peak response ---
    [~, peakLag] = max(max(max(abs(strf))));
    peakFrame = strf(:,:,peakLag);
    [~, peakPix] = max(abs(peakFrame(:)));
    [py, px] = ind2sub([40 40], peakPix);
    
    x0_init   = px;
    y0_init   = py;
    sig_init  = 5;
    amp_init  = squeeze(mean(mean(strf, 1), 2))';  % (1 x 32) mean amplitude per lag
    
    params0 = [x0_init, y0_init, sig_init, amp_init];
    
    % --- Bounds ---
    lb = [1,  1,  1,  -inf*ones(1,32)];
    ub = [40, 40, 20,  inf*ones(1,32)];
    
    % --- Fit ---
    opts = optimset('Display', 'off', 'MaxIter', 500);
    params_fit = lsqcurvefit(gaussModel, params0, xdata, ydata, lb, ub, opts);
    
    allParams(n,:) = params_fit;
    
    % Reconstruct fitted STRF
    fitted_vec = gaussModel(params_fit, xdata);
    allFitted{n} = reshape(fitted_vec, 40, 40, 32);
end

%% Part B — Plot empirical vs fitted STRF at 80ms + amplitude time course

lags = 300:-10:-10;
lagIdx80 = find(lags == 80);
figure;

for n = 1:length(neurons)
    neuron = neurons(n);
    strf = plot_strf(retinaData, neuron);
    fitted = allFitted{n};
    amplitudes = squeeze(allParams(n, 4:end));  % 32 amplitudes
    
    % --- Subplot 1: empirical STRF at 80ms ---
    subplot(length(neurons), 3, (n-1)*3 + 1);
    max_val = max(abs(strf(:)));
    img = strf(:,:,lagIdx80);
    image((img + max_val) / (2*max_val) * 255);
    axis image off;
    if n == 1, title('Empirical (80ms)'); end
    ylabel(sprintf('N%d', neuron), 'FontSize', 8);
    
    % --- Subplot 2: fitted model at 80ms ---
    subplot(length(neurons), 3, (n-1)*3 + 2);
    img_fit = fitted(:,:,lagIdx80);
    max_fit = max(abs(fitted(:)));
    image((img_fit + max_fit) / (2*max_fit) * 255);
    axis image off;
    if n == 1, title('Fitted model (80ms)'); end
    
    % --- Subplot 3: amplitude time course ---
    subplot(length(neurons), 3, (n-1)*3 + 3);
    plot(lags, amplitudes, 'k', 'LineWidth', 1.2);
    set(gca, 'XDir', 'reverse');  % 300ms on left → 0ms on right
    xlabel('Lag (ms)'); 
    yline(0, '--', 'Color', [0.6 0.6 0.6]);
    xlim([-10 300]);
    if n == 1, title('Amplitude time course'); end
    box off;
end
colormap(gray(256));
sgtitle('STRF Gaussian Fits — All Neurons');

%Part C
%% Part C — Summary statistics

sigmas       = allParams(:, 3);          % receptive field widths
amplitudes   = allParams(:, 4:end);      % (10 x 32) amplitude time courses

% Average receptive field width (in pixels — convert if pixel size is known)
avg_width = mean(sigmas);
fprintf('Average RF width (sigma): %.2f pixels\n', avg_width);

% OFF cells: peak amplitude is negative (activated by darkening)
[~, peakLagIdx] = max(abs(amplitudes), [], 2);   % index of peak lag per neuron
peakAmps = arrayfun(@(i) amplitudes(i, peakLagIdx(i)), 1:length(neurons))';

isOFF = peakAmps < 0;
pctOFF = 100 * sum(isOFF) / length(neurons);
fprintf('Percentage OFF cells: %.1f%%\n', pctOFF);

% Average time to peak for OFF and ON cells
peakTimes = lags(peakLagIdx)';   % ms

avg_peak_OFF = mean(peakTimes(isOFF));
avg_peak_ON  = mean(peakTimes(~isOFF));
fprintf('Average time to peak — OFF cells: %.1f ms\n', avg_peak_OFF);
fprintf('Average time to peak — ON  cells: %.1f ms\n', avg_peak_ON);