function pcaTrajectories(R)
% pcaTrajectories  Run PCA trajectory analysis on neural data.
%
%   pcaTrajectories(R)  takes the trial struct array R and produces:
%     1. A figure showing the top 3 PCs over time for straight vs. curved
%        reaches under three normalization schemes (none, z-score, soft).
%     2. A 3D figure of neural trajectories grouped by reach direction.
%
%   R must contain fields:
%     R(t).conditionID       – condition identifier
%     R(t).moveOnsetTime     – movement-onset time (ms)
%     R(t).unit(u).spikeTimes – spike times for unit u (ms)
%     R(t).targetXY          – 2-element [x; y] target position

% -------------------------------------------------------------------------
% Parameters
% -------------------------------------------------------------------------
pre_window  = 500;
post_window = 1500;
bin_size    = 20;
k           = 5;       % number of CV folds
soft_norm   = 5;       % soft normalization constant (for CV section)
soft_norm_alpha = 5;   % soft normalization constant (for direction PCA)

methods            = {'none', 'zscore', 'soft'};
normalization_labels = {'No Normalization', 'Z-Score', 'Soft Normalization'};

% -------------------------------------------------------------------------
% Trial indexing
% -------------------------------------------------------------------------
conditionIDs   = [R.conditionID];
straight_idx   = mod(conditionIDs, 3) == 1;
curved_idx     = mod(conditionIDs, 3) == 2;
straight_trials = find(straight_idx);
curved_trials   = find(curved_idx);

time_bins   = -pre_window : bin_size : post_window;
bin_centers = time_bins(1:end-1) + bin_size/2;
n_bins      = numel(bin_centers);
n_units     = numel(R(1).unit);

onset_align = [R.moveOnsetTime];
all_trials  = [straight_trials(:); curved_trials(:)];
n_trials    = numel(all_trials);

% -------------------------------------------------------------------------
% Build firing rate tensor  [n_trials x n_units x n_bins]
% -------------------------------------------------------------------------
firing_rates = nan(n_trials, n_units, n_bins);
for t = 1:n_trials
    tr = all_trials(t);
    t0 = onset_align(tr);
    for u = 1:n_units
        spike_times = R(tr).unit(u).spikeTimes - t0;
        counts = histcounts(spike_times, time_bins);
        firing_rates(t, u, :) = counts / (bin_size / 1000);
    end
end
disp('Firing rates created.')

% -------------------------------------------------------------------------
% Condition-averaged PSTHs  [2 x n_bins x n_units]
% -------------------------------------------------------------------------
fr_straight  = squeeze(mean(firing_rates(1:numel(straight_trials), :, :), 1));
fr_curved    = squeeze(mean(firing_rates(numel(straight_trials)+1:end, :, :), 1));
fr_cond_avg  = permute(cat(3, fr_straight', fr_curved'), [3 1 2]);
n_conds      = 2;

% -------------------------------------------------------------------------
% Cross-validated PCA dimensionality analysis
% -------------------------------------------------------------------------
n_straight   = numel(straight_trials);
n_curved     = numel(curved_trials);
cond_labels  = [ones(n_straight, 1); 2*ones(n_curved, 1)];

fold_ids_trials = zeros(n_trials, 1);
for cond = 1:2
    cond_mask = find(cond_labels == cond);
    n_cond    = numel(cond_mask);
    perm      = randperm(n_cond);
    cond_folds = mod(0:n_cond-1, k) + 1;
    fold_ids_trials(cond_mask) = cond_folds(perm);
end

n_X_fixed  = floor(n_units / 2);
max_dims   = n_X_fixed - 1;
n_plot_dims = min(n_X_fixed - 1, 2*n_bins - 1);

cv_loss      = nan(numel(methods), k, k, max_dims);
all_explained = nan(numel(methods), n_plot_dims);
all_score    = cell(numel(methods), 1);
all_coeff    = cell(numel(methods), 1);

for m = 1:numel(methods)

    % --- Cross-validation loop ---
    for fold_r = 1:k
        train_idx = find(fold_ids_trials ~= fold_r);
        test_idx  = find(fold_ids_trials == fold_r);
        n_train   = numel(train_idx);
        n_test    = numel(test_idx);

        for fold_c = 1:k
            neuron_perm = randperm(n_units);
            x_neurons   = neuron_perm(1:n_X_fixed);
            y_neurons   = neuron_perm(n_X_fixed+1:end);
            n_Y         = numel(y_neurons);

            X_train = reshape(permute(firing_rates(train_idx, x_neurons, :), [1 3 2]), n_train*n_bins, n_X_fixed);
            Y_train = reshape(permute(firing_rates(train_idx, y_neurons, :), [1 3 2]), n_train*n_bins, n_Y);
            X_test  = reshape(permute(firing_rates(test_idx,  x_neurons, :), [1 3 2]), n_test*n_bins,  n_X_fixed);
            Y_test  = reshape(permute(firing_rates(test_idx,  y_neurons, :), [1 3 2]), n_test*n_bins,  n_Y);

            switch methods{m}
                case 'none'
                    X_train_n = X_train; X_test_n = X_test;
                case 'zscore'
                    mu_x  = mean(X_train, 1);
                    sig_x = std(X_train, 0, 1); sig_x(sig_x == 0) = 1;
                    X_train_n = (X_train - mu_x) ./ sig_x;
                    X_test_n  = (X_test  - mu_x) ./ sig_x;
                case 'soft'
                    rng_x     = max(X_train, [], 1) - min(X_train, [], 1);
                    X_train_n = X_train ./ (rng_x + soft_norm);
                    X_test_n  = X_test  ./ (rng_x + soft_norm);
            end

            mu_X      = mean(X_train_n, 1);
            X_train_c = X_train_n - mu_X;
            X_test_c  = X_test_n  - mu_X;
            mu_Y      = mean(Y_train, 1);
            Y_train_c = Y_train - mu_Y;

            [coeff, score, ~] = pca(X_train_c, 'NumComponents', max_dims);

            for d = 1:max_dims
                S_train = score(:, 1:d);
                beta    = S_train \ Y_train_c;
                S_test  = X_test_c * coeff(:, 1:d);
                Y_hat   = S_test * beta + mu_Y;
                cv_loss(m, fold_r, fold_c, d) = mean((Y_test - Y_hat).^2, 'all');
            end
        end
    end

    % --- Full-data PCA for visualization ---
    FR_all = reshape(permute(firing_rates, [1 3 2]), n_trials*n_bins, n_units);

    switch methods{m}
        case 'none'
            FR_all_n = FR_all;
        case 'zscore'
            mu_x  = mean(FR_all, 1);
            sig_x = std(FR_all, 0, 1); sig_x(sig_x == 0) = 1;
            FR_all_n = (FR_all - mu_x) ./ sig_x;
        case 'soft'
            rng_x    = max(FR_all, [], 1) - min(FR_all, [], 1);
            FR_all_n = FR_all ./ (rng_x + soft_norm);
    end

    FR_all_n  = reshape(FR_all_n, n_trials, n_bins, n_units);
    psth_s    = squeeze(mean(FR_all_n(1:n_straight, :, :), 1));
    psth_c    = squeeze(mean(FR_all_n(n_straight+1:end, :, :), 1));
    X_plot_final = [psth_s; psth_c];
    X_plot_final = X_plot_final - mean(X_plot_final, 1);

    [coeff_full, score_full, ~, ~, explained] = pca(X_plot_final);
    all_score{m}        = score_full;
    all_explained(m, :) = explained(1:n_plot_dims)';
    all_coeff{m}        = coeff_full;
end

mean_loss = squeeze(mean(cv_loss, [2 3])); %#ok<NASGU>

% =========================================================================
% Figure 1 – Top 3 PCs over time (straight vs. curved, 3 normalizations)
% =========================================================================
figure('Position', [100 100 1300 500]);
pc_colors = nebula(3);   % replace with nebula(3) if that colormap is available

for m = 1:numel(methods)
    subplot(1, 3, m); hold on;
    score = all_score{m};
    for pc = 1:3
        plot(bin_centers, score(1:n_bins,    pc), '-',  'Color', pc_colors(pc,:), 'LineWidth', 2);
        plot(bin_centers, score(n_bins+1:end, pc), '--', 'Color', pc_colors(pc,:), 'LineWidth', 2);
    end
    xline(0, 'k--', 'LineWidth', 1.5);
    xlabel('Time from movement onset (ms)');
    ylabel('PC score');
    title(normalization_labels{m});
    if m == 1
        legend('PC1 straight','PC1 curved', ...
               'PC2 straight','PC2 curved', ...
               'PC3 straight','PC3 curved', ...
               'Location','best');
    end
end
sgtitle('Top 3 PCs Over Time');

% =========================================================================
% Figure 2 – 3D neural trajectories by reach direction
% =========================================================================
targetXY    = reshape([R.targetXY], 2, [])';
angles      = atan2d(targetXY(:,2), targetXY(:,1));
angleGroups = round(angles / 45) * 45;
angleGroups(angleGroups == 180)              = -180;
angleGroups(angleGroups == 0 & angles <  5)  = -45;
angleGroups(angleGroups == 0 & angles >= 5)  =  45;
uniqueGroups = unique(angleGroups);
nGroups      = length(uniqueGroups);
dirLabels    = {'Left','Lower-Left','Lower-Right','Upper-Right','Up','Upper-Left'};

all_trial_ids = 1:numel(R);
n_all         = numel(all_trial_ids);

fr_all = nan(n_all, n_units, n_bins);
for t = 1:n_all
    tr = all_trial_ids(t);
    t0 = onset_align(tr);
    for u = 1:n_units
        spike_times        = R(tr).unit(u).spikeTimes - t0;
        counts             = histcounts(spike_times, time_bins);
        fr_all(t, u, :)   = counts / (bin_size / 1000);
    end
end

psth_groups = nan(n_units, nGroups * n_bins);
for g = 1:nGroups
    group_trials = find(angleGroups == uniqueGroups(g));
    psth_groups(:, (g-1)*n_bins + (1:n_bins)) = ...
        squeeze(mean(fr_all(group_trials, :, :), 1));
end

fr_range = max(psth_groups, [], 2) - min(psth_groups, [], 2);
X_dir    = psth_groups ./ (fr_range + soft_norm_alpha);
X_dir    = X_dir - mean(X_dir, 2);

[~, score_dir, ~, ~, explained_dir] = pca(X_dir'); %#ok<ASGLU>

dir_colors = hsv(nGroups);
[~, onset_bin] = min(abs(bin_centers));

figure('Name', '3D Neural Trajectories', 'Position', [100 100 900 700]);
hold on;
h_legend = gobjects(nGroups, 1);

for g = 1:nGroups
    idx = (g-1)*n_bins + (1:n_bins);
    pc1 = score_dir(idx, 1);
    pc2 = score_dir(idx, 2);
    pc3 = score_dir(idx, 3);
    col = dir_colors(g, :);

    n_seg  = n_bins - 1;
    alphas = linspace(0.25, 1.0, n_seg);
    for s = 1:n_seg
        seg_col = col .* alphas(s) + (1 - alphas(s)) * [1 1 1];
        seg_col = min(max(seg_col, 0), 1);
        plot3(pc1(s:s+1), pc2(s:s+1), pc3(s:s+1), '-', ...
            'Color', seg_col, 'LineWidth', 2.5, 'HandleVisibility', 'off');
    end

    h_legend(g) = plot3(nan, nan, nan, '-', ...
        'Color', col, 'LineWidth', 2, ...
        'DisplayName', sprintf('%d°', uniqueGroups(g)));
    plot3(pc1(1),         pc2(1),         pc3(1),         '^', ...
        'MarkerFaceColor', col, 'MarkerEdgeColor', 'k', ...
        'MarkerSize', 7,  'HandleVisibility', 'off');
    plot3(pc1(onset_bin), pc2(onset_bin), pc3(onset_bin), 'o', ...
        'MarkerFaceColor', col, 'MarkerEdgeColor', 'k', ...
        'MarkerSize', 8,  'HandleVisibility', 'off');
end

plot3(0, 0, 0, 'k+', 'MarkerSize', 12, 'LineWidth', 2, 'DisplayName', 'Origin');
xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
title('3D neural trajectories by reach direction (soft norm, all-data PCA)');
legend(dirLabels, 'Location', 'best', 'NumColumns', 2);
annotation('textbox', [0.01 0.01 0.3 0.06], ...
    'String', '▲ = start   ● = movement onset', ...
    'EdgeColor', 'none', 'FontSize', 9, 'Color', [0.4 0.4 0.4]);
grid on; axis equal; view(35, 25); box off;
set(gca, 'FontSize', 12);

end