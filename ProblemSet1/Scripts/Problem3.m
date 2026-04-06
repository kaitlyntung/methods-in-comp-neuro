%% Load data and set variables
load('Data\retinaData.mat');

%% 3.1 Spike Triggered Average
% Part A
function strf = plot_strf(retinaData, neuron)
    % data cleaning, get rid of all that have 0 like problem set says to do
    valid = retinaData.stimulusFrameTimes ~= 0;
    stimulusFrameTimes = retinaData.stimulusFrameTimes(valid);
    stimulus = retinaData.stimulusFrames(:,:,valid) - 0.5;
    stimulus = double(stimulus);
    stimulusTimes = stimulusFrameTimes * 0.1; % convert to ms
    stimulusTimes = double(stimulusTimes);
    spikeTimes = retinaData.spikes{neuron} * 0.1;
    spikeTimes = double(spikeTimes);
    
    lag_times = 300:-10:-10;
    numLags = length(lag_times);
    numSpikes = length(spikeTimes);
    
    stimulus_reshaped = reshape(stimulus, 1600, []);
    spike_idx = round(interp1(stimulusTimes, ...
        1:length(stimulusTimes), ...
        spikeTimes, ...
        'nearest', ...
        'extrap'));
    spike_idx = max(1, min(spike_idx, length(stimulusTimes)));
    
    strf = zeros(40, 40, numLags);
    for lag = 1:numLags
        lagFrames = round(lag_times(lag) / (1000/30));
        offset_idx = spike_idx - lagFrames;
        valid_idx = offset_idx >= 1 & offset_idx <= size(stimulus_reshaped, 2);
        strf(:,:,lag) = reshape(sum(stimulus_reshaped(:, offset_idx(valid_idx)), 2), 40, 40);
    end
    strf = strf / numSpikes;
    
    % Plotting
    figure;
    max_val = max(abs(strf(:)));
    for k = 1:numLags
        subplot(6, 6, k);
        img_scaled = (strf(:,:,k) + max_val) / (2 * max_val);
        image(img_scaled * 255);
        axis image off;
        title(sprintf('%d ms', lag_times(k)), 'FontSize', 7);
    end
    colormap(gray(256));
    sgtitle(sprintf('Neuron %d', neuron));
end

%Part B
neurons = [1, 4, 11, 15, 26, 51, 80, 84, 96, 105];
strfs = cell(length(neurons), 1);

for i = 1:length(neurons)
    strfs{i} = plot_strf(retinaData, neurons(i));
end

%% 3.2 Receptive Field Models

neurons = [1, 4, 11, 15, 26, 51, 80, 84, 96, 105];
lag_times = 300:-10:-10; 
numLags = length(lag_times);
lagIdx80 = find(lag_times == 80);

allParams  = zeros(length(neurons), 3 + numLags);  % [x0, y0, sigma, A_1...A_32]
allFitted  = cell(length(neurons), 1);
allSTRFs   = cell(length(neurons), 1);

%% Part A – Fit time-varying 2D Gaussian with lsqcurvefit
% Model: shared (x0, y0, sigma), one amplitude per lag
% params = [x0, y0, sigma, A_1, A_2, ..., A_32]  (length = 3 + 32 = 35)
%
% xdata passed to lsqcurvefit must be a matrix; ydata must be a vector.
% Strategy: replicate the 40x40 pixel grid for each lag, stack into one
% long vector of length 1600*32 = 51200.

[X, Y] = meshgrid(1:40, 1:40);
xyGrid  = [X(:), Y(:)];                    % 1600 x 2, one row per pixel
% Replicate grid across all lags so lsqcurvefit sees one big xdata matrix
xdata_full = repmat(xyGrid, numLags, 1);   % (1600*32) x 2

gaussModel = @(params, xdata) ...
    cell2mat( ...
        arrayfun(@(k) ...
            params(3+k) * exp( ...
                -( (xdata(1:1600,1)-params(1)).^2 + ...
                   (xdata(1:1600,2)-params(2)).^2 ) ...
                / (2*params(3)^2) ), ...
        1:numLags, 'UniformOutput', false)' ...
    );
% Note: because every 1600-row block of xdata_full is identical (same pixel
% grid), we can simply always index rows 1:1600 for the spatial part and
% vary only the amplitude params(3+k).  This avoids indexing into xdata
% inside the loop which would require passing lag boundaries separately.

for n = 1:length(neurons)
    neuron = neurons(n);
    fprintf('Fitting neuron %d (%d/%d)...\n', neuron, n, length(neurons));

    strf = plot_strf(retinaData, neuron);   % your existing function
    close(gcf);                             % suppress the figure it opens
    allSTRFs{n} = strf;

    % ydata: stack all lag frames into one long column vector
    % strf is (40 x 40 x 32); reshape to (1600 x 32) then to (51200 x 1)
    ydata = reshape(strf, 1600, numLags);   % 1600 x 32
    ydata = ydata(:);                       % 51200 x 1  (pixel-major, lag-minor)

    % ---- Warm start ------------------------------------------------
    [~, peakLag] = max(squeeze(max(max(abs(strf)))));  % lag with largest response
    peakFrame    = strf(:,:,peakLag);
    [~, peakPix] = max(abs(peakFrame(:)));
    [py, px]     = ind2sub([40 40], peakPix);

    x0_init  = double(px);
    y0_init  = double(py);
    sig_init = 5;
    % One amplitude per lag: use the spatial mean of each lag frame
    amp_init = squeeze(mean(mean(strf, 1), 2))';   % 1 x 32

    params0 = [x0_init, y0_init, sig_init, amp_init];   % 1 x 35

    % ---- Bounds ----------------------------------------------------
    lb = [1,  1,  0.5, -inf*ones(1,numLags)];
    ub = [40, 40, 20,   inf*ones(1,numLags)];

    % ---- Fit -------------------------------------------------------
    opts = optimset('Display', 'off', 'MaxIter', 800, 'TolFun', 1e-8);

    params_fit = lsqcurvefit(gaussModel, params0, xdata_full, ydata, lb, ub, opts);

    allParams(n,:) = params_fit;

    % ---- Reconstruct fitted STRF -----------------------------------
    fitted_vec  = gaussModel(params_fit, xdata_full);           % 51200 x 1
    allFitted{n} = reshape(fitted_vec, 40, 40, numLags);
end

%% Part B – Plot: empirical STRF | model STRF | amplitude time course
figure('Position', [50 50 1400 900]);
for n = 1:length(neurons)
    neuron   = neurons(n);
    strf     = allSTRFs{n};
    fitted   = allFitted{n};
    amps     = allParams(n, 4:end);   % A_1 … A_32

    % ---- Shared color scale (empirical frame at 80 ms) ----
    empFrame    = strf(:,:,lagIdx80);
    fitFrame    = fitted(:,:,lagIdx80);
    clim_val    = max(abs([empFrame(:); fitFrame(:)]));
    if clim_val == 0, clim_val = 1; end

    % ---- Subplot 1: empirical STRF at 80 ms -------------------
    subplot(length(neurons), 3, (n-1)*3 + 1);
    imagesc(empFrame, [-clim_val, clim_val]);
    axis image off;
    colormap(gca, 'gray');
    if n == 1, title('Empirical (80 ms)', 'FontSize', 8); end
    ylabel(sprintf('N%d', neuron), 'FontSize', 7);

    % ---- Subplot 2: fitted STRF at 80 ms ----------------------
    subplot(length(neurons), 3, (n-1)*3 + 2);
    imagesc(fitFrame, [-clim_val, clim_val]);
    axis image off;
    colormap(gca, 'gray');
    if n == 1, title('Model (80 ms)', 'FontSize', 8); end

    % ---- Subplot 3: amplitude time course ---------------------
    subplot(length(neurons), 3, (n-1)*3 + 3);
    plot(lag_times, amps, 'k-', 'LineWidth', 1.2);
    xlim([min(lag_times) max(lag_times)]);
    set(gca, 'XDir', 'reverse');   % 300 ms on the left → 0 on the right
    yline(0, 'r--', 'LineWidth', 0.8);
    xlabel('Lag (ms)', 'FontSize', 6);
    ylabel('Amplitude', 'FontSize', 6);
    if n == 1, title('Amp. time course', 'FontSize', 8); end
end
sgtitle('3.2B – Empirical vs. Fitted STRF and Amplitude Time Courses');

%% Part C – Summary statistics

sigma_vals   = allParams(:, 3);               % Gaussian width for each neuron
avgWidth     = mean(sigma_vals);

% Identify ON vs OFF cells:
%   OFF cell → peak amplitude is NEGATIVE (responds to darkening)
%   ON  cell → peak amplitude is POSITIVE (responds to brightening)
peakAmps = zeros(length(neurons), 1);
peakTimes = zeros(length(neurons), 1);
for n = 1:length(neurons)
    amps = allParams(n, 4:end);               % 1 x 32
    [~, pkIdx] = max(abs(amps));              % index of largest |amplitude|
    peakAmps(n)  = amps(pkIdx);
    peakTimes(n) = lag_times(pkIdx);          % corresponding lag in ms
end

isOFF = peakAmps < 0;
isON  = peakAmps > 0;

pctOFF       = 100 * sum(isOFF) / length(neurons);
avgPeakOFF   = mean(peakTimes(isOFF));
avgPeakON    = mean(peakTimes(isON));

fprintf('\n===== 3.2C Summary =====\n');
fprintf('Average receptive field width (sigma): %.2f pixels\n', avgWidth);
fprintf('OFF cells: %d / %d  (%.1f%%)\n', sum(isOFF), length(neurons), pctOFF);
fprintf('Average time-to-peak – OFF cells: %.1f ms\n', avgPeakOFF);
fprintf('Average time-to-peak – ON  cells: %.1f ms\n', avgPeakON);