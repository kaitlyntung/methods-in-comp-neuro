% We will compute the Fano Factor as a function of time aligned to both the 
% target appearance, movement onset, and go-cue separately for curved and straight 
% reaches to determine how variability evolves across delay and movement 
% and whether this differs for the reach types, considering curved reaches 
% are arguably more complex.

conditionIDs = [R.conditionID];
straight_idx = mod(conditionIDs, 3) == 1;
curved_idx = mod(conditionIDs, 3) == 2;

straight_trials = find(straight_idx);
curved_trials   = find(curved_idx);

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

straight_dir_ids = floor(conditionIDs(straight_idx) / 3); 
curved_dir_ids = floor(conditionIDs(curved_idx) / 3);

unique_straight_dirs = unique(straight_dir_ids);
unique_curved_dirs = unique(curved_dir_ids);

ff_target = nan(num_units, num_bins, 2);
ff_gocue = nan(num_units, num_bins, 2);
ff_onset = nan(num_units, num_bins, 2);

for unit = 1:num_units
    for cond = 1:2
        if cond == 1
            trial_ids = straight_trials;
            dir_ids = straight_dir_ids;
            unique_dirs = unique_straight_dirs;
        else
            trial_ids = curved_trials;
            dir_ids = curved_dir_ids;
            unique_dirs = unique_curved_dirs;
        end

        num_dirs = numel(unique_dirs);

        ff_target_by_dir = nan(num_dirs, num_bins);
        ff_gocue_by_dir = nan(num_dirs, num_bins);
        ff_onset_by_dir = nan(num_dirs, num_bins);

        for d = 1:num_dirs
            dir_mask = dir_ids == unique_dirs(d);
            dir_trial_ids = trial_ids(dir_mask); 
            num_dir_trials = numel(dir_trial_ids);

            counts_target = nan(num_dir_trials, num_bins);
            counts_gocue  = nan(num_dir_trials, num_bins);
            counts_onset  = nan(num_dir_trials, num_bins);

            for trial = 1:num_dir_trials
                tr = dir_trial_ids(trial);
                spike_times = R(tr).unit(unit).spikeTimes;

                rel_target = spike_times - target_align(tr);
                rel_gocue = spike_times - gocue_align(tr);
                rel_onset = spike_times - onset_align(tr);

                for b = 1:num_bins
                    t_lo = bin_centers(b) - bin_size/2;
                    t_hi = bin_centers(b) + bin_size/2;

                    counts_target(trial, b) = sum(rel_target >= t_lo & rel_target < t_hi);
                    counts_gocue( trial, b) = sum(rel_gocue  >= t_lo & rel_gocue  < t_hi);
                    counts_onset( trial, b) = sum(rel_onset  >= t_lo & rel_onset  < t_hi);
                end
            end

            ff_target_by_dir(d, :) = compute_fano_factor(counts_target);
            ff_gocue_by_dir( d, :) = compute_fano_factor(counts_gocue);
            ff_onset_by_dir( d, :) = compute_fano_factor(counts_onset);
        end

        ff_target(unit, :, cond) = nanmean(ff_target_by_dir, 1);
        ff_gocue( unit, :, cond) = nanmean(ff_gocue_by_dir,  1);
        ff_onset( unit, :, cond) = nanmean(ff_onset_by_dir,  1);
    end
end

mean_ff = @(ff, cond) nanmean(ff(:, :, cond), 1);
sem_ff  = @(ff, cond) nanstd(ff(:, :, cond), 0, 1) / sqrt(sum(~all(isnan(ff(:,:,cond)), 2)));

%Plotting
colors = struct('straight', [0 0 0.8], 'curved', [0.8 0 0]);
align_labels = {'Target Appearance', 'Go Cue', 'Movement Onset'};
x_labels     = {'Time from Target Appearance (ms)', ...
                 'Time from Go Cue (ms)', ...
                 'Time from Movement Onset (ms)'};
ff_all = {ff_target, ff_gocue, ff_onset};

figure;

for a = 1:3
    subplot(1, 3, a);
    hold on;

    ff = ff_all{a};

    for cond = 1:2
        m  = mean_ff(ff, cond);
        sem = sem_ff(ff, cond);
        if cond == 1
            col = colors.straight;
            lbl = 'Straight';
            fill_lbl = 'SEM Straight';
        else
            col = colors.curved;
            lbl = 'Curved';
            fill_lbl = 'SEM Curved';
        end

        fill([bin_centers, fliplr(bin_centers)], ...
             [m + sem, fliplr(m - sem)], ...
             col, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', fill_lbl);

        plot(bin_centers, m, '-', 'Color', col, 'LineWidth', 2, 'DisplayName', lbl);
    end

    xline(0, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    yline(1, 'k:',  'LineWidth', 1.0, 'HandleVisibility', 'off');

    xlabel(x_labels{a});
    ylabel('Fano Factor');
    title(align_labels{a});
    xlim([-pre_window, post_window]);
    if a == 1
        legend('Location', 'best');
    end
end

sgtitle('Population Fano Factor');

%% Helper function
function ff = compute_fano_factor(counts)
    m  = mean(counts, 1);
    v  = var(counts,  0, 1);
    ff = nan(size(m));
    valid = m > 0;
    ff(valid) = v(valid) ./ m(valid);
end
