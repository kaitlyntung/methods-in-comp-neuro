% -------------------------------------------------------------------------
% Fano Factor as a function of time: aligned to target appearance and
% movement onset, split by curved vs straight reaches
% -------------------------------------------------------------------------

% Time windows
bin_size    = 50;   % ms
pre_window  = 500;  % ms before alignment event
post_window = 800;  % ms after alignment event
time_bins   = -pre_window : bin_size : post_window;
bin_centers = time_bins(1:end-1) + bin_size/2;

% Identify valid trial indices for each condition
straight_trials = find(straight_idx);
curved_trials   = find(curved_idx);

% -------------------------------------------------------------------------
% Helper: count spikes in bins for a single unit across trials
% aligned to a given event time per trial
% -------------------------------------------------------------------------
function counts = bin_spikes(spike_times, align_times, trial_ids, time_bins, R)
    n_trials = numel(trial_ids);
    n_bins   = numel(time_bins) - 1;
    counts   = nan(n_trials, n_bins);
    for t = 1:n_trials
        tr        = trial_ids(t);
        t0        = align_times(tr);
        rel_times = spike_times - t0;
        for b = 1:n_bins
            counts(t, b) = sum(rel_times >= time_bins(b) & ...
                               rel_times <  time_bins(b+1));
        end
    end
end

% -------------------------------------------------------------------------
% Compute Fano Factor: var / mean across trials at each time bin
% Fano is undefined (set to NaN) when mean == 0
% -------------------------------------------------------------------------
function ff = compute_fano(counts)
    m  = mean(counts, 1);
    v  = var(counts,  0, 1);
    ff = nan(size(m));
    ff(m > 0) = v(m > 0) ./ m(m > 0);
end

% -------------------------------------------------------------------------
% Loop over units
% -------------------------------------------------------------------------
n_units = numel(R(1).unit);

% Preallocate: units x time bins x 2 (straight/curved) x 2 (target/onset)
ff_target = nan(n_units, numel(bin_centers), 2);
ff_onset  = nan(n_units, numel(bin_centers), 2);

target_align = [R.targetAppearsTime];
onset_align  = [R.moveOnsetTime];

for u = 1:n_units

    % Collect spike times for this unit across all trials
    % Assumes R(trial).unit(u).spikeTimes holds spike times in ms
    for cond = 1:2
        if cond == 1
            trial_ids = straight_trials;
        else
            trial_ids = curved_trials;
        end

        n_trials = numel(trial_ids);
        n_bins   = numel(bin_centers);

        counts_target = nan(n_trials, n_bins);
        counts_onset  = nan(n_trials, n_bins);

        for t = 1:n_trials
            tr          = trial_ids(t);
            spike_times = R(tr).unit(u).spikeTimes;

            % Align to target appearance
            t0_target = target_align(tr);
            rel_target = spike_times - t0_target;

            % Align to movement onset
            t0_onset  = onset_align(tr);
            rel_onset  = spike_times - t0_onset;

            for b = 1:n_bins
                counts_target(t, b) = sum(rel_target >= time_bins(b) & ...
                                          rel_target <  time_bins(b+1));
                counts_onset(t, b)  = sum(rel_onset  >= time_bins(b) & ...
                                          rel_onset  <  time_bins(b+1));
            end
        end

        ff_target(u, :, cond) = compute_fano(counts_target);
        ff_onset(u,  :, cond) = compute_fano(counts_onset);
    end
end

% -------------------------------------------------------------------------
% Average Fano Factor across units
% -------------------------------------------------------------------------
mean_ff_target_straight = nanmean(ff_target(:, :, 1), 1);
mean_ff_target_curved   = nanmean(ff_target(:, :, 2), 1);
mean_ff_onset_straight  = nanmean(ff_onset(:,  :, 1), 1);
mean_ff_onset_curved    = nanmean(ff_onset(:,  :, 2), 1);

sem_ff_target_straight  = nanstd(ff_target(:, :, 1), 0, 1) / sqrt(n_units);
sem_ff_target_curved    = nanstd(ff_target(:, :, 2), 0, 1) / sqrt(n_units);
sem_ff_onset_straight   = nanstd(ff_onset(:,  :, 1), 0, 1) / sqrt(n_units);
sem_ff_onset_curved     = nanstd(ff_onset(:,  :, 2), 0, 1) / sqrt(n_units);

% -------------------------------------------------------------------------
% Plot
% -------------------------------------------------------------------------
colors = [0.29 0.47 0.81;   % straight (blue)
          0.84 0.37 0.37];  % curved   (red)

figure;
set(gcf, 'Color', 'w', 'Position', [100 100 1100 450]);

align_labels = {'Target Appearance', 'Movement Onset'};

for panel = 1:2
    subplot(1, 2, panel);
    hold on;

    if panel == 1
        ff_straight = mean_ff_target_straight;
        ff_curved   = mean_ff_target_curved;
        se_straight = sem_ff_target_straight;
        se_curved   = sem_ff_target_curved;
    else
        ff_straight = mean_ff_onset_straight;
        ff_curved   = mean_ff_onset_curved;
        se_straight = sem_ff_onset_straight;
        se_curved   = sem_ff_onset_curved;
    end

    % Shaded SEM
    fill([bin_centers, fliplr(bin_centers)], ...
         [ff_straight + se_straight, fliplr(ff_straight - se_straight)], ...
         colors(1,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    fill([bin_centers, fliplr(bin_centers)], ...
         [ff_curved + se_curved, fliplr(ff_curved - se_curved)], ...
         colors(2,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none');

    % Mean lines
    plot(bin_centers, ff_straight, 'Color', colors(1,:), 'LineWidth', 2);
    plot(bin_centers, ff_curved,   'Color', colors(2,:), 'LineWidth', 2);

    % Alignment marker
    xline(0, 'k--', 'LineWidth', 1.5);

    xlabel(sprintf('Time relative to %s (ms)', align_labels{panel}), 'FontSize', 11);
    ylabel('Fano Factor', 'FontSize', 11);
    title(sprintf('Fano Factor aligned to\n%s', align_labels{panel}), ...
          'FontSize', 12, 'FontWeight', 'bold');
    legend({'', 'No Barriers', '', 'Barriers'}, 'Location', 'best');
    box off;
end

sgtitle('Trial-to-Trial Variability: Fano Factor by Condition', ...
        'FontSize', 13, 'FontWeight', 'bold');