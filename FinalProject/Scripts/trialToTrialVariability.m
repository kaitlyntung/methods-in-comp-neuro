% We will compute the Fano Factor as a function of time aligned to both the 
% target appearance and movement onset separately for curved and straight 
% reaches to determine how variability evolves across delay and movement 
% and whether this differs for the reach types, considering curved reaches 
% are arguably more complex.
conditionIDs = [R.conditionID];
straight_idx = mod(conditionIDs, 3) == 1;
curved_idx = mod(conditionIDs, 3) == 2;

straight_trials = find(straight_idx);
curved_trials = find(curved_idx);

bin_size = 20;
pre_window = 200;
post_window = 750;
time_bins = -pre_window : bin_size : post_window;
bin_centers = time_bins(1:end-1) + bin_size/2;

num_units = numel(R(1).unit); % channels

ff_target = nan(num_units, numel(bin_centers), 2);
ff_onset = nan(num_units, numel(bin_centers), 2);

target_align = [R.targetAppearsTime];
onset_align = [R.moveOnsetTime];

for unit = 1:num_units
    for cond = 1:2
        if cond == 1
            trial_ids = straight_trials;
        else
            trial_ids = curved_trials;
        end
        num_trials = numel(trial_ids);
        num_bins = numel(bin_centers);

        counts_target = nan(num_trials, num_bins);
        counts_onset = nan(num_trials, num_bins);

        for trial = 1:num_trials
            tr = trial_ids(trial);
            spike_times = R(tr).unit(unit).spikeTimes;

            time_target = target_align(tr);
            time_onset = onset_align(tr);
            relative_to_target = spike_times - time_target;
            relative_to_onset = spike_times - time_onset;

            for bins = 1:num_bins
                counts_target(trial, bins) = sum(relative_to_target >= time_bins(bins) & ...
                                          relative_to_target <  time_bins(bins+1));
                counts_onset(trial, bins) = sum(relative_to_onset  >= time_bins(bins) & ...
                                          relative_to_onset  <  time_bins(bins+1));
            end
        end

        ff_target(unit, :, cond) = compute_fano_factor(counts_target);
        ff_onset(unit,  :, cond) = compute_fano_factor(counts_onset);
    end
end

mean_ff_target_straight = nanmean(ff_target(:, :, 1), 1);
mean_ff_target_curved = nanmean(ff_target(:, :, 2), 1);
mean_ff_onset_straight = nanmean(ff_onset(:,  :, 1), 1);
mean_ff_onset_curved = nanmean(ff_onset(:,  :, 2), 1);

sem_ff_target_straight = nanstd(ff_target(:, :, 1), 0, 1) / sqrt(num_units);
sem_ff_target_curved = nanstd(ff_target(:, :, 2), 0, 1) / sqrt(num_units);
sem_ff_onset_straight = nanstd(ff_onset(:,  :, 1), 0, 1) / sqrt(num_units);
sem_ff_onset_curved = nanstd(ff_onset(:,  :, 2), 0, 1) / sqrt(num_units);

%% Plotting
figure;

% Target-aligned
subplot(1,2,1);
hold on;

fill([bin_centers, fliplr(bin_centers)], ...
     [mean_ff_target_straight + sem_ff_target_straight, fliplr(mean_ff_target_straight - sem_ff_target_straight)], ...
     [0 0 0.8], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
fill([bin_centers, fliplr(bin_centers)], ...
     [mean_ff_target_curved + sem_ff_target_curved, fliplr(mean_ff_target_curved - sem_ff_target_curved)], ...
     'red', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

plot(bin_centers, mean_ff_target_straight, 'b-', 'LineWidth', 2);
plot(bin_centers, mean_ff_target_curved,   'r-', 'LineWidth', 2);
xline(0, 'k--', 'LineWidth', 1.5);
yline(1, 'k:', 'LineWidth', 1);

xlabel('Time from Target Appearance (ms)');
ylabel('Fano Factor');
title('Target-Aligned');
legend('No Barriers ± SEM', 'Barriers ± SEM', 'Location', 'best');

% Onset-aligned
subplot(1,2,2);
hold on;

fill([bin_centers, fliplr(bin_centers)], ...
     [mean_ff_onset_straight + sem_ff_onset_straight, fliplr(mean_ff_onset_straight - sem_ff_onset_straight)], ...
     [0 0 0.8], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
fill([bin_centers, fliplr(bin_centers)], ...
     [mean_ff_onset_curved + sem_ff_onset_curved, fliplr(mean_ff_onset_curved - sem_ff_onset_curved)], ...
     'red', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

plot(bin_centers, mean_ff_onset_straight, 'b-', 'LineWidth', 2);
plot(bin_centers, mean_ff_onset_curved,   'r-', 'LineWidth', 2);

xline(0, 'k--', 'LineWidth', 1.5);
yline(1, 'k:', 'LineWidth', 1);

xlabel('Time from Movement Onset (ms)');
ylabel('Fano Factor');
title('Movement-Aligned');
legend('No Barriers', 'Barriers', 'Location', 'best');

sgtitle('Trial by Trial Variability (Fano Factor)');

%% Helper functions
function ff = compute_fano_factor(counts)
    m = mean(counts, 1);
    v = var(counts,  0, 1);
    ff = nan(size(m));
    ff(m > 0) = v(m > 0) ./ m(m > 0);
end