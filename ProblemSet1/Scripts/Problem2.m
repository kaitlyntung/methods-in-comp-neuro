%% Load data and set variables
load('Data\mtSpikeTimes.mat');

%% 2.1 Basic Spike Train Statistics
% Part A
num_trials = length(mtSpikeTimes);
count_30 = zeros(num_trials, 1);

for i = 1:num_trials
    spike_times = mtSpikeTimes{i} * 1000; % convert to ms
    if isempty(spike_times)
        continue;
    end
    first_spike = spike_times(1);
    spikes_after = spike_times(spike_times > first_spike);
    count_30(i) = sum(spikes_after <= first_spike + 30);
end
firing_rate_30 = count_30/0.03;
avg_firing_rate_30 = mean(firing_rate_30);

fprintf('Average firing rate (30-ms window): %d \n', avg_firing_rate_30);

% Part B
fano_factor_30 = var(counts_30) / mean(counts_30);
fprintf('Fano factor (30-ms window): %d \n', fano_factor_30);

% Part C
count_300 = zeros(num_trials, 1);

for i = 1:num_trials
    spike_times = mtSpikeTimes{i} * 1000; % convert to ms
    if isempty(spike_times)
        continue;
    end
    first_spike = spike_times(1);
    spikes_after = spike_times(spike_times > first_spike);
    count_300(i) = sum(spikes_after <= first_spike + 30);
end
firing_rate_300 = count_300/0.3;
avg_firing_rate_300 = mean(firing_rate_300);

fprintf('Average firing rate (300-ms window): %d \n', avg_firing_rate_300);
fano_factor_300 = var(counts_300) / mean(counts_300);
fprintf('Fano factor (300-ms window): %d \n', fano_factor_300);

%% 2.2 Interspike Intervals
all_isi = [];

for i = 1:length(mtSpikeTimes)
    spike_times = mtSpikeTimes{i} * 1000; % ms
    if isempty(spike_times)
        continue;
    end
    first_spike = spike_times(1);
    spike_range = spike_times(spike_times > first_spike & spike_times <= first_spike + 300);

    if length(spike_range) > 1
        isi = diff(spike_range);
        all_isi = [all_isi, isi];
    end
end

figure;
hold on;
bins = 0:10:300;
counts = histcounts(all_isi, bins);
bin_centers = bins(1:end-1) + diff(bins)/2;
bar(bin_centers, counts, 'FaceColor', 'red');

% Possion calculations
lambda = avg_firing_rate_300;
t = bin_centers / 1000; % convert from ms to seconds
pdf = lambda * exp(-lambda * t);

bin_width = 10 / 1000; % 10ms bins like we did above
expected_vals = pdf * length(all_isi) * bin_width;
plot(bin_centers, expected_vals, 'k', 'LineWidth', 2);

set(gca, 'FontName', 'Arial');    
xlabel('Time (ms)');
ylabel('Count');
title('Interspike Interval Histogram and Expected Value');

legend('ISI', 'Expected Value (Poisson)');
box off;