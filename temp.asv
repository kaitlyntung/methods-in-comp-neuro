%% Setup and firing rate computation
conditionIDs = [R.conditionID];
straight_idx = mod(conditionIDs, 3) == 1;
curved_idx   = mod(conditionIDs, 3) == 2;

straight_trials = find(straight_idx);
curved_trials   = find(curved_idx);

pre_window  = 500;
post_window = 1500;
bin_size    = 20;
time_bins   = -pre_window : bin_size : post_window;
bin_centers = time_bins(1:end-1) + bin_size/2;
n_bins      = numel(bin_centers);
n_units     = numel(R(1).unit);

onset_align = [R.moveOnsetTime];
all_trials  = [straight_trials(:); curved_trials(:)];
n_trials    = numel(all_trials);

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
disp('firing rates created')

n_straight = numel(straight_trials);
n_curved   = numel(curved_trials);

psth_straight = squeeze(mean(firing_rates(1:n_straight, :, :), 1));
psth_curved   = squeeze(mean(firing_rates(n_straight+1:end, :, :), 1));
X_psth = [psth_straight, psth_curved]'; 

col_range  = max(X_psth, [], 1) - min(X_psth, [], 1);
soft_denom = col_range + 5;

mu_unit    = mean(X_psth, 1);
sigma_unit = std(X_psth, 0, 1);
sigma_unit(sigma_unit == 0) = 1;

normalization_labels = {'No Normalization', 'Z-Score', 'Soft Normalization'};
methods       = {'none', 'zscore', 'soft'};
colors_method = [0 0 0; 0 0 0.8; 0.8 0 0];
soft_norm_alpha = 5;
n_folds  = 5;
max_dims = 20;
n_cols   = size(X_psth, 1);

cv_loss     = nan(numel(methods), max_dims);
all_coeff   = cell(numel(methods), 1);
all_score   = cell(numel(methods), 1);
all_explained = nan(numel(methods), n_units);

for m = 1:numel(methods)

    switch methods{m}
        case 'none'
            X_norm = X_psth;
        case 'zscore'
            X_norm = (X_psth - mu_unit) ./ sigma_unit;
        case 'soft'
            X_norm = X_psth ./ soft_denom;
    end

    X_norm = X_norm - mean(X_norm, 1);

    [coeff, score, ~, ~, explained] = pca(X_norm);
    all_coeff{m}   = coeff;
    all_score{m}   = score;
    all_explained(m, 1:numel(explained)) = explained;

    idx       = randperm(n_cols);
    cv_idx    = zeros(n_cols, 1);
    fold_size = floor(n_cols / n_folds);
    for f = 1:n_folds
        cv_idx(idx((f-1)*fold_size + 1 : f*fold_size)) = f;
    end
    cv_idx(cv_idx == 0) = n_folds;

    fold_loss = nan(n_folds, max_dims);
    for f = 1:n_folds
        test_mask  = (cv_idx == f);
        train_mask = ~test_mask;

        X_train = X_norm(train_mask, :);
        X_test  = X_norm(test_mask,  :); 

        train_mean = mean(X_train, 1);
        X_train    = X_train - train_mean;
        X_test     = X_test  - train_mean;

        W_all = pca(X_train);

        for k = 1:max_dims
            W     = W_all(:, 1:k);
            X_recon = (X_test * W) * W';
            fold_loss(f, k) = mean((X_test(:) - X_recon(:)).^2);
        end
    end

    cv_loss(m, :) = mean(fold_loss, 1);
end

%% loss curve
figure; hold on;
for m = 1:numel(methods)
    loss      = cv_loss(m, :);
    loss_norm = (loss - min(loss)) / (max(loss) - min(loss));
    plot(1:max_dims, loss_norm, '-o', ...
        'Color', colors_method(m,:), 'LineWidth', 2);
end
xlabel('Number of PCs (k)');
ylabel('Normalized Reconstruction Error (0-1)');
title('Cross-Validated PCA Loss by Dimensionality');
legend(normalization_labels, 'Location', 'northeast');
grid on;

%% scree plot
n_plot = 20;
figure; hold on;
for m = 1:numel(methods)
    plot(1:n_plot, all_explained(m, 1:n_plot), '-o', ...
        'Color', colors_method(m,:), 'LineWidth', 2);
end
xlabel('Principal Component');
ylabel('Variance Explained (%)');
title('Scree Plot');
legend(normalization_labels, 'Location', 'northeast');
grid on;

%% cumulative variance
figure; hold on;
for m = 1:numel(methods)
    cum_exp = cumsum(all_explained(m, :));
    plot(1:n_plot, cum_exp(1:n_plot), '-o', ...
        'Color', colors_method(m,:), 'LineWidth', 2);
end
yline(80, 'k--', 'LineWidth', 1.5, 'Label', '80%');
yline(95, 'k:',  'LineWidth', 1.5, 'Label', '95%');
xlabel('Number of Components');
ylabel('Cumulative Variance Explained (%)');
title('Cumulative Variance');
ylim([0 100]);
legend(normalization_labels, 'Location', 'southeast');
grid on;

%% 2d pc trajectories (over time)
figure;
pc_colors = nebula(3);
for m = 1:numel(methods)
    subplot(1, 3, m);
    hold on;
    score = all_score{m};
    for k = 1:3
        pc_straight = score(1:n_bins, k);
        pc_curved   = score(n_bins+1:end, k);
        plot(bin_centers, pc_straight, '-',  'Color', pc_colors(k,:), 'LineWidth', 2);
        plot(bin_centers, pc_curved,   '--', 'Color', pc_colors(k,:), 'LineWidth', 2);
    end
    xline(0, 'k--', 'LineWidth', 1.5);
    xlabel('Time from Movement Onset (ms)');
    ylabel('PC Score');
    title(normalization_labels{m});
    if m == 1
        legend('PC1 Straight', 'PC1 Curved', 'PC2 Straight', 'PC2 Curved', ...
               'PC3 Straight', 'PC3 Curved', 'Location', 'best');
    end
end
sgtitle('Top 3 PCs Over Time by Normalization Method');

%% 3d trajectories (by direction) 
targetXY    = reshape([R.targetXY], 2, [])';
angles      = atan2d(targetXY(:,2), targetXY(:,1));
angleGroups = round(angles / 45) * 45;
angleGroups(angleGroups == 180) = -180;
uniqueGroups = unique(angleGroups);
nGroups      = length(uniqueGroups);

all_trial_ids = 1:numel(R);
n_all         = numel(all_trial_ids);

fr_all = nan(n_all, n_units, n_bins);
for t = 1:n_all
    tr = all_trial_ids(t);
    t0 = onset_align(tr);
    for u = 1:n_units
        spike_times = R(tr).unit(u).spikeTimes - t0;
        counts = histcounts(spike_times, time_bins);
        fr_all(t, u, :) = counts / (bin_size / 1000);
    end
end

% Build [n_units x nGroups*n_bins] PSTH matrix
psth_groups = nan(n_units, nGroups * n_bins);
for g = 1:nGroups
    group_mask   = angleGroups == uniqueGroups(g);
    group_trials = find(group_mask);
    psth_groups(:, (g-1)*n_bins + (1:n_bins)) = ...
        squeeze(mean(fr_all(group_trials, :, :), 1));
end

% Z-score normalization for trajectory plot (transpose to [obs x units] first)
X_dir  = psth_groups';  % [nGroups*n_bins x n_units]
mu_dir    = mean(X_dir, 1);
sig_dir   = std(X_dir, 0, 1);
sig_dir(sig_dir == 0) = 1;
X_dir_norm = (X_dir - mu_dir) ./ sig_dir;
X_dir_norm = X_dir_norm - mean(X_dir_norm, 1);

[~, score_dir, ~, ~, explained_dir] = pca(X_dir_norm);

dir_colors = hsv(nGroups);
[~, onset_bin] = min(abs(bin_centers));

figure('Position', [100, 100, 900, 750]);
hold on;

h_legend = gobjects(nGroups, 1);
for g = 1:nGroups
    idx = (g-1)*n_bins + (1:n_bins);
    pc1 = score_dir(idx, 1);
    pc2 = score_dir(idx, 2);
    pc3 = score_dir(idx, 3);
    col = dir_colors(g, :);

    h_legend(g) = plot3(pc1, pc2, pc3, '-', ...
        'Color', col, 'LineWidth', 2, ...
        'DisplayName', sprintf('%d°', uniqueGroups(g)));

    % Movement onset marker
    plot3(pc1(onset_bin), pc2(onset_bin), pc3(onset_bin), 'o', ...
        'MarkerFaceColor', col, 'MarkerEdgeColor', 'k', ...
        'MarkerSize', 8, 'HandleVisibility', 'off');

    % Start marker
    plot3(pc1(1), pc2(1), pc3(1), '^', ...
        'MarkerFaceColor', col, 'MarkerEdgeColor', 'k', ...
        'MarkerSize', 7, 'HandleVisibility', 'off');

    % Fading color gradient over time
    n_seg  = n_bins - 1;
    alphas = linspace(0.25, 1.0, n_seg);
    for s = 1:n_seg
        seg_col = col .* alphas(s) + (1 - alphas(s)) * [1 1 1];
        seg_col = min(max(seg_col, 0), 1);
        plot3(pc1(s:s+1), pc2(s:s+1), pc3(s:s+1), '-', ...
            'Color', seg_col, 'LineWidth', 2.5, 'HandleVisibility', 'off');
    end
end

plot3(0, 0, 0, 'k+', 'MarkerSize', 12, 'LineWidth', 2, 'DisplayName', 'Origin');

xlabel(sprintf('PC1 (%.1f%%)', explained_dir(1)));
ylabel(sprintf('PC2 (%.1f%%)', explained_dir(2)));
zlabel(sprintf('PC3 (%.1f%%)', explained_dir(3)));
title('3D Neural Trajectories by Reach Direction (Z-scored)');
legend(h_legend, 'Location', 'bestoutside', 'NumColumns', 2);
annotation('textbox', [0.01, 0.01, 0.3, 0.06], ...
    'String', '▲ = start   ● = movement onset', ...
    'EdgeColor', 'none', 'FontSize', 9, 'Color', [0.4 0.4 0.4]);
grid on; axis equal; view(35, 25);