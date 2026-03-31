%% Load data and set variables
load('Data\mtSpikeTimes.mat');

%% 1.1 Raster Plots
figure;
hold on

num_trials = length(mtSpikeTimes);

spike_height = 0.8; 

for trial = 1:num_trials
    if ~isempty(mtSpikeTimes{trial})
        spike_times = mtSpikeTimes{trial} * 1000; % convert to ms

        for s = 1:length(spike_times)
            x = spike_times(s);
            line([x x], ...
                 [trial - spike_height/2, trial + spike_height/2], ...
                 'Color', 'k', 'LineWidth', 1.5);
        end
    end
end

set(gca, 'FontName', 'Arial');
title('\bf Raster Plot', 'FontSize', 18)

text(0, 75, '\bf Trial Number', 'FontSize', 12, 'Rotation', 90)
text(225, -7, '\bf Time since trial start (ms)', 'FontSize', 12)

ylim([-5, num_trials])
xlim([0, 750])

line([0, 0], [0, 20], 'color', 'k', 'LineWidth', 4)
text(0, 20, '\bf 20', 'FontSize', 10, 'Rotation', 90)

line([0, 50], [0, 0], 'color', 'k', 'LineWidth', 2)
text(50, 0, '\bf 50', 'FontSize', 10)

axis off
box off
%% 1.2 Peri-Stimulus Time Histograms
% Part A
function func = psth(data, bin_size)
    bins = 0:bin_size:750;
    all_spike_times = [data{:}]; % gets all trials concatenated
    all_spike_times = all_spike_times * 1000;
    spike_counts = histcounts(all_spike_times, bins);
    num_trials = length(data);
    firingRate = spike_counts / ((bin_size / 1000) * num_trials);
    figure;
    b = bar(bins(1:end-1), firingRate, 'histc');
    b.FaceColor = 'red';
    
    set(gca, 'FontName', 'Arial');
    title(sprintf('PSTH Pooled Across Trials (%dms bins)', bin_size), 'FontSize', 18)
    text(0, 20, '\bf Firing Rate (spikes/s)', 'FontSize', 12, 'Rotation', 90)
    text(300, -2, '\bf Time since trial start (ms)', 'FontSize', 12)

    line([0, 0], [10, 0], 'color', 'k', 'LineWidth', 4)
    text(0, 10, '\bf 10', 'FontSize', 10, 'Rotation', 90)
    line([0, 100], [0, 0], 'color', 'k', 'LineWidth', 4)
    text(100, -2, '\bf 100', 'FontSize', 10)
    axis off
    box off
end

psth(mtSpikeTimes, 10);
psth(mtSpikeTimes, 20);

% Part B
function psth_line(data, bin_size)

    bins = 0:bin_size:750;
    num_trials = length(data);
    num_bins = length(bins) - 1;
    trial_rates = zeros(num_trials, num_bins);

    for i = 1:num_trials
        spike_times = data{i} * 1000; % convert to ms
        spike_counts = histcounts(spike_times, bins);
        trial_rates(i, :) = spike_counts / (bin_size / 1000);
    end

    mean_rate = mean(trial_rates, 1);
    sem_rate = std(trial_rates, 0, 1) / sqrt(num_trials);
    bin_centers = bins(1:end-1) + bin_size/2;

    figure;
    hold on;

    fill([bin_centers fliplr(bin_centers)], ...
         [mean_rate + sem_rate, fliplr(mean_rate - sem_rate)], ...
         'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    plot(bin_centers, mean_rate, 'r', 'LineWidth', 2);

    title(sprintf('PSTH Line Plot (%d ms bins)', bin_size), 'FontSize', 18)
    xlabel('Time since trial start (ms)')
    ylabel('Firing Rate (spikes/s)')

    set(gca, 'FontName', 'Arial');    
    text(0, 20, '\bf Firing Rate (spikes/s)', 'FontSize', 12, 'Rotation', 90)
    text(300, -2, '\bf Time since trial start (ms)', 'FontSize', 12)

    line([0, 0], [10, 0], 'color', 'k', 'LineWidth', 4)
    text(0, 10, '\bf 10', 'FontSize', 10, 'Rotation', 90)
    line([0, 100], [0, 0], 'color', 'k', 'LineWidth', 4)
    text(100, -2, '\bf 100', 'FontSize', 10)
    axis off
    box off
end
psth_line(mtSpikeTimes, 10);
psth_line(mtSpikeTimes, 20);