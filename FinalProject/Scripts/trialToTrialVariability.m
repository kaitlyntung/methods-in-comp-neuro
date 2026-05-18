% We will compute the Fano Factor as a function of time aligned to both the 
% target appearance, movement onset, and go-cue separately for curved and straight 
% reaches to determine how variability evolves across delay and movement 
% and whether this differs for the reach types, considering curved reaches 
% are arguably more complex.

conditionIDs = [R.conditionID];
straight_idx = mod(conditionIDs, 3) == 1;
curved_idx = mod(conditionIDs, 3) == 2;

straight_trials = find(straight_idx);
curved_trials = find(curved_idx);

bin_size = 50;
step_size = 10;
pre_window = 400;
post_window = 600;

bin_centers = (-pre_window + bin_size/2) : step_size : (post_window - bin_size/2);
num_bins = numel(bin_centers);
num_units = numel(R(1).unit);

target_align = [R.targetAppearsTime];
gocue_align = [R.goCueTime];
onset_align = [R.moveOnsetTime];

ff_target = nan(num_units, num_bins, 2);
ff_gocue = nan(num_units, num_bins, 2);
ff_onset = nan(num_units, num_bins, 2);

for unit = 1:num_units
    for cond = 1:2
        if cond == 1
            trial_ids = straight_trials;
        else
            trial_ids = curved_trials;
        end
        num_trials = numel(trial_ids);

        counts_target = nan(num_trials, num_bins);
        counts_gocue = nan(num_trials, num_bins);
        counts_onset = nan(num_trials, num_bins);

        for trial = 1:num_trials
            tr = trial_ids(trial);
            spike_times = R(tr).unit(unit).spikeTimes;

            rel_target = spike_times - target_align(tr);
            rel_gocue = spike_times - gocue_align(tr);
            rel_onset = spike_times - onset_align(tr);

            for b = 1:num_bins
                t_lo = bin_centers(b) - bin_size/2;
                t_hi = bin_centers(b) + bin_size/2;

                counts_target(trial, b) = sum(rel_target >= t_lo & rel_target < t_hi);
                counts_gocue(trial, b) = sum(rel_gocue >= t_lo & rel_gocue < t_hi);
                counts_onset(trial, b) = sum(rel_onset >= t_lo & rel_onset < t_hi);
            end
        end

        ff_target(unit, :, cond) = compute_fano_factor(counts_target);
        ff_gocue(unit, :, cond) = compute_fano_factor(counts_gocue);
        ff_onset(unit, :, cond) = compute_fano_factor(counts_onset);
    end
end

mean_ff = @(ff, cond) nanmean(ff(:, :, cond), 1);
sem_ff = @(ff, cond) nanstd(ff(:, :, cond), 0, 1) / sqrt(sum(~all(isnan(ff(:,:,cond)), 2)));

%% Plotting
colors = struct('straight', [0 0 0.8], 'curved', [0.8 0 0]);
align_labels = {'Target Appearance', 'Go Cue', 'Movement Onset'};
x_labels     = {'Time from Target Appearance (ms)', ...
                 'Time from Go Cue (ms)', ...
                 'Time from Movement Onset (ms)'};
ff_all = {ff_target, ff_gocue, ff_onset};

figure('Position', [100, 100, 1400, 450]);

for a = 1:3
    subplot(1, 3, a);
    hold on;

    ff = ff_all{a};

    for cond = 1:2
        m   = mean_ff(ff, cond);
        sem = sem_ff(ff, cond);

        if cond == 1
            col  = colors.straight;
            lbl  = 'No Barriers';
        else
            col  = colors.curved;
            lbl  = 'Barriers';
        end

        % Shaded SEM band
        fill([bin_centers, fliplr(bin_centers)], ...
             [m + sem, fliplr(m - sem)], ...
             col, 'FaceAlpha', 0.2, 'EdgeColor', 'none');

        % Mean line
        plot(bin_centers, m, '-', 'Color', col, 'LineWidth', 2, 'DisplayName', lbl);
    end

    xline(0, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    yline(1, 'k:',  'LineWidth', 1.0, 'HandleVisibility', 'off');

    xlabel(x_labels{a});
    ylabel('Fano Factor');
    title(align_labels{a});
    % legend('No Barriers', 'Barriers', 'Location', 'best');
    xlim([-pre_window, post_window]);
end

sgtitle('Population Fano Factor');

%% Helper function
function ff = compute_fano_factor(counts)
    % counts: trials x bins
    % Returns FF only for bins where mean > 0; NaN otherwise
    m  = mean(counts, 1);
    v  = var(counts,  0, 1);
    ff = nan(size(m));
    valid     = m > 0;
    ff(valid) = v(valid) ./ m(valid);
end