%% Problem 1
numTimes = size(mtNeuron.data, 1);
numDirections = size(mtNeuron.data, 2);
numTrials = size(mtNeuron.data, 3);

colors = hsv(numDirections);
directions = linspace(-90, 90, numDirections);

figure; hold on;
trial_offset = 0;
legend_handles = gobjects(numDirections,1);

for direction = 1:numDirections
    for trial = 1:numTrials
        % convert to ms
        spikes = find(mtNeuron.data(:, direction, trial)) * 2;
        y = (trial_offset + trial) * ones(size(spikes));

        plot(spikes, y, '.', 'Color', colors(direction, :), 'MarkerSize', 1);
    end
    
    legend_handles(direction) = plot(nan, nan, '.', ...
        'Color', colors(direction, :), 'MarkerSize', 10);
    
    trial_offset = trial_offset + numTrials;
end

xlabel('Time (ms)');
ylabel('Trial Number');
title('Raster Plot');
legend(legend_handles, string(directions) + "°", 'Location', 'eastoutside');
%% Problem 2
cumulative_spikes = cumsum(mtNeuron.data, 1);
num_directions = size(mtNeuron.data, 2);
num_trials = size(mtNeuron.data, 3);
num_time = size(mtNeuron.data, 1);
% I(T) = I(T_n, theta) = sum_n P_T(n)log_s(P_T(n|theta)/P_T(n))

fractions = [1/4, 1/3, 1/2, 2/3, 1]; % eq7
num_fracs = length(fractions);
MI_raw = zeros(num_time, num_fracs);

for fi = 1:num_fracs
    n_samp = floor(fractions(fi) * num_trials); 
    for t = 1:num_time

        counts = zeros(num_directions, n_samp);
        for d = 1:num_directions
            idx = randperm(num_trials, n_samp);
            counts(d, :) = cumulative_spikes(t, d, idx);
        end

        all_counts = counts(:);
        max_n = max(all_counts);
        count_vals = 0:max_n;
        num_vals = length(count_vals);

        % P(n|theta) based on histogram
        Pn_given_theta = zeros(num_directions, num_vals);
        for d = 1:num_directions
            for ni = 1:num_vals
                Pn_given_theta(d, ni) = mean(counts(d, :) == count_vals(ni));
            end
        end
        % P(n) = (1/num_directions) * sum_theta P(n|theta) eq2, avg
        Pn = mean(Pn_given_theta, 1);
        MI = 0;
        P_theta = 1 / num_directions; 

        for d = 1:num_directions
            % P(n | theta_d)
            Pn_d = Pn_given_theta(d, :);    
            mask = Pn_d > 0 & Pn > 0;
            MI = MI + P_theta * sum( Pn_d(mask) .* log2( Pn_d(mask) ./ Pn(mask)) );
        end

        MI_raw(t, fi) = MI;
    end
end
MI_corrected = zeros(num_time, 1);
inv_fracs    = 1 ./ fractions;

for t = 1:num_time
    p = polyfit(inv_fracs, MI_raw(t, :), 1);  
    MI_corrected(t) = p(2); % y-intercept = I_true
end

MI_corrected = max(MI_corrected, 0);
time_ms = (0:num_time-1) * 2;

figure;
plot(time_ms, MI_corrected, 'k-', 'LineWidth', 1.5);
xlabel('Time (ms)');
ylabel('MI (bits)');
title('Mutual Information');
xlim([0 512]);

%% Problem 3
n_boot = 1000;                        
MI_boot = zeros(num_time, n_boot);

for b = 1:n_boot
    boot_spikes = zeros(num_time, num_directions, num_trials);
    for d = 1:num_directions
        idx = randi(num_trials, 1, num_trials);
        boot_spikes(:, d, :) = cumulative_spikes(:, d, idx);
    end

    MI_raw_boot = zeros(num_time, num_fracs);

    for fi = 1:num_fracs
        n_samp = floor(fractions(fi) * num_trials);

        for t = 1:num_time
            counts = zeros(num_directions, n_samp);
            for d = 1:num_directions
                idx = randperm(num_trials, n_samp);
                counts(d, :) = boot_spikes(t, d, idx);
            end

            all_counts = counts(:);
            max_n = max(all_counts);
            count_vals = 0:max_n;
            num_vals = length(count_vals);

            Pn_given_theta = zeros(num_directions, num_vals);
            for d = 1:num_directions
                for ni = 1:num_vals
                    Pn_given_theta(d, ni) = mean(counts(d, :) == count_vals(ni));
                end
            end

            Pn = mean(Pn_given_theta, 1);
            MI = 0;
            P_theta = 1 / num_directions;
            for d = 1:num_directions
                Pn_d = Pn_given_theta(d, :);
                mask = Pn_d > 0 & Pn > 0;
                MI = MI + P_theta * sum(Pn_d(mask) .* log2(Pn_d(mask) ./ Pn(mask)));
            end
            MI_raw_boot(t, fi) = MI;
        end
    end

    % Bias correction for this bootstrap replicate
    for t = 1:num_time
        p = polyfit(inv_fracs, MI_raw_boot(t, :), 1);
        MI_boot(t, b) = p(2);
    end
    MI_boot(:, b) = max(MI_boot(:, b), 0);
end

% Percentile-based confidence interval
MI_low  = prctile(MI_boot, 2.5,  2);   % 2.5th percentile across boots
MI_high = prctile(MI_boot, 97.5, 2);   % 97.5th percentile

%% Plotting
figure;
fill([time_ms, fliplr(time_ms)], ...
     [MI_low', fliplr(MI_high')], ...
     [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
hold on;
plot(time_ms, MI_corrected, 'k-', 'LineWidth', 1.5);
xlabel('Time (ms)');
ylabel('MI (bits)');
title('Mutual Information with Bootstrap');
xlim([0 512]);
legend('95% Bootstrapped Confidence Interval', 'Mutual Information', 'Location', 'best');

%% Problem 4 — Latency detection and MI fractions

% ── Method: baseline threshold ──────────────────────────────────────────────
% We define latency as the first time bin at which the bias-corrected MI
% exceeds the pre-stimulus mean + 2 SD of the pre-stimulus MI.
% The pre-stimulus period is t < 0, i.e. bins before stimulus onset (t=0).
% Using the data's own noise floor is principled: it adapts to this neuron's
% variability rather than relying on an arbitrary absolute threshold.
%
% The stimulus begins at t=0 (bin 1). Since recordings start at t=0 and
% the array has no negative-time baseline, we use the first ~50ms (bins 1-25)
% as a proxy pre-stimulus period — before the neuron could plausibly respond
% (given typical MT latencies of ~50-80ms). We require the threshold to be
% exceeded for at least 3 consecutive bins to avoid single-bin noise spikes.

% ── Define proxy pre-stimulus window ────────────────────────────────────────
pre_bins    = 1:25;          % 0–50 ms: before any plausible neural response
baseline_mu = mean(MI_corrected(pre_bins));
baseline_sd = std(MI_corrected(pre_bins));
threshold   = baseline_mu + 2 * baseline_sd;

% ── Find first sustained crossing ───────────────────────────────────────────
n_consec = 3;   % require 3 consecutive bins above threshold (~6 ms)
above    = MI_corrected > threshold;
latency_bin = NaN;
for t = 1:(num_time - n_consec + 1)
    if all(above(t:t+n_consec-1))
        latency_bin = t;
        break;
    end
end
latency_ms = (latency_bin - 1) * 2;   % convert bin → ms
fprintf('Estimated latency: %d ms (bin %d)\n', latency_ms, latency_bin);

% ── MI fractions ─────────────────────────────────────────────────────────────
% "Within the first X ms of the neural response" means X ms after latency,
% not after stimulus onset.

total_MI = MI_corrected(end);   % final (asymptotic) value

% 50 ms after latency
bin_50  = latency_bin + round(50  / 2);   % 50ms / 2ms per bin
bin_100 = latency_bin + round(100 / 2);   % 100ms / 2ms per bin
bin_50  = min(bin_50,  num_time);
bin_100 = min(bin_100, num_time);

MI_at_50  = MI_corrected(bin_50);
MI_at_100 = MI_corrected(bin_100);

frac_50  = MI_at_50  / total_MI;
frac_100 = MI_at_100 / total_MI;

fprintf('MI at 50ms after latency:  %.3f bits (%.1f%% of total)\n', MI_at_50,  frac_50*100);
fprintf('MI at 100ms after latency: %.3f bits (%.1f%% of total)\n', MI_at_100, frac_100*100);

% ── Plot ─────────────────────────────────────────────────────────────────────
figure;
fill([time_ms, fliplr(time_ms)], ...
     [MI_low', fliplr(MI_high')], ...
     [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
hold on;
plot(time_ms, MI_corrected, 'k-', 'LineWidth', 1.5);
yline(threshold,    'b--', 'Threshold',          'LabelVerticalAlignment','bottom');
xline(latency_ms,   'r--', sprintf('Latency: %d ms', latency_ms));
xline(time_ms(bin_50),  'g--', '+50ms');
xline(time_ms(bin_100), 'm--', '+100ms');

% Annotate fractions
text(time_ms(bin_50)+4,  MI_at_50,  sprintf('%.1f%%', frac_50*100),  'Color','g','FontSize',9);
text(time_ms(bin_100)+4, MI_at_100, sprintf('%.1f%%', frac_100*100), 'Color','m','FontSize',9);

xlabel('Time (ms)');
ylabel('MI (bits)');
title('Mutual Information with Bootstrap');
xlim([0 512]);
legend('95% Bootstrapped Confidence Interval', 'Mutual Information', 'Location', 'best');
