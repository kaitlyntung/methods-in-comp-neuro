%% Part 1
%% Part A
nConditions = length(R);
nUnits  = length(R(1).unit);
startTime = -150;
endTime = 500;
binSize = 10;
bin_edges = (startTime:binSize:endTime);
nTimeBins = length(bin_edges) - 1; 
spikeCounts = zeros(nTimeBins, nConditions, nUnits);

for trial = 1:nConditions
    moveOnsetTime = R(trial).moveOnsetTime;
    for unit = 1:nUnits
        spikeTime = R(trial).unit(unit).spikeTimes;
        if isempty(spikeTime)
            continue
        end
        spikeTimeLocked = spikeTime - moveOnsetTime;
        spikeCounts(:, trial, unit) = histcounts(spikeTimeLocked, bin_edges);
    end
end

%% Part B
sd_ms = 40;
sd_bins = sd_ms / binSize;
spike_counts_smooth = gaussFilt1(spikeCounts, sd_bins, 1);

%% Part C
% Convert from spikes/bin to spikes/s at the same time.
% spikes/s = spikes/bin * (1000 ms/s / bin_size ms)

conditionIDs = [R.conditionID];
conditions   = unique(conditionIDs);   % sorted list of condition IDs
nConditions  = length(conditions);

% Output: nTimeBins x nConditions x nUnits
mean_fr = zeros(nTimeBins, nConditions, nUnits);

for c = 1:nConditions
    trial_idx = conditionIDs == conditions(c);   % logical index into trials
    % Mean over matching trials (dim 2), then convert to spikes/s
    mean_fr(:, c, :) = mean(spike_counts_smooth(:, trial_idx, :), 2) ...
                       * (1000 / binSize);
end

fprintf('mean_fr size: %d x %d x %d  (spikes/s)\n', size(mean_fr));

%% 6. Plot PSTHs for units 1, 4, 101, 104

units_to_plot = [1, 4, 101, 104];
nPlotUnits    = length(units_to_plot);

figure('Name', 'PSTHs', 'Position', [100 100 1000 700]);

for p = 1:nPlotUnits
    u = units_to_plot(p);

    subplot(2, 2, p);
    hold on;

    % Plot one trace per condition; color by condition index
    cmap = lines(nConditions);
    for c = 1:nConditions
        plot(binCenters, mean_fr(:, c, u), 'Color', cmap(c, :), 'LineWidth', 1);
    end

    xline(0, 'k--', 'LineWidth', 1.5);   % mark movement onset
    xlabel('Time from movement onset (ms)');
    ylabel('Firing rate (spikes/s)');
    title(sprintf('Unit %d', u));
    xlim([bin_centers(1) bin_centers(end)]);
    box off;
end

sgtitle('PSTHs locked to movement onset');