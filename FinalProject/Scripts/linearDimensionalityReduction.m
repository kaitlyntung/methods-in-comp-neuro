% -------------------------------------------------------------------------
% Cross-validated PCA on M1 activity
% Comparing dimensionality under no normalization, z-score, soft normalization
% -------------------------------------------------------------------------

% Time window for neural activity (relative to movement onset, ms)
pre_window  = 500;
post_window = 800;
bin_size    = 50;
time_bins   = -pre_window : bin_size : post_window;
bin_centers = time_bins(1:end-1) + bin_size/2;
n_bins      = numel(bin_centers);
n_units     = numel(R(1).unit);

onset_align = [R.moveOnsetTime];
all_trials  = [straight_trials(:); curved_trials(:)];
n_trials    = numel(all_trials);

% -------------------------------------------------------------------------
% Step 1: Build firing rate matrix [trials x units x time bins]
% -------------------------------------------------------------------------
FR = nan(n_trials, n_units, n_bins);

for t = 1:n_trials
    tr = all_trials(t);
    for u = 1:n_units
        spike_times = R(tr).unit(u).spikeTimes;
        t0          = onset_align(tr);
        rel_times   = spike_times - t0;
        for b = 1:n_bins
            FR(t, u, b) = sum(rel_times >= time_bins(b) & ...
                               rel_times <  time_bins(b+1)) ...
                           / (bin_size / 1000);  % convert to Hz
        end
    end
end

% -------------------------------------------------------------------------
% Step 2: Define normalization schemes
% -------------------------------------------------------------------------
norm_labels = {'No Normalization', 'Z-Score', 'Soft Normalization'};
soft_norm_alpha = 5;  % spikes/s, controls soft normalization floor

normalize_FR = @(X, method) apply_normalization(X, method, soft_norm_alpha);

function X_norm = apply_normalization(X, method, alpha)
    % X is [trials x units x bins]; normalize across trials per unit per bin
    % Reshape to [trials x (units*bins)] for convenience
    [n_trials, n_units, n_bins] = size(X);
    X_flat = reshape(X, n_trials, n_units * n_bins);

    switch method
        case 'none'
            X_norm = X_flat;

        case 'zscore'
            mu     = mean(X_flat, 1);
            sigma  = std(X_flat,  0, 1);
            sigma(sigma == 0) = 1;  % avoid divide by zero for silent units
            X_norm = (X_flat - mu) ./ sigma;

        case 'soft'
            % Soft normalization: divide by (range + alpha)
            fr_range = max(X_flat, [], 1) - min(X_flat, [], 1);
            X_norm   = X_flat ./ (fr_range + alpha);
    end

    X_norm = reshape(X_norm, n_trials, n_units, n_bins);
end

% -------------------------------------------------------------------------
% Step 3: Cross-validated PCA
% Randomly split trials into two halves. Fit PCA on half A, project half B.
% Compute cross-validated variance explained as covariance between
% projections of the two halves onto the PCs from half A.
% Repeat n_cv times and average.
% -------------------------------------------------------------------------
n_cv    = 50;    % number of cross-validation folds
n_comps = min(n_units, n_trials) - 1;

cv_var_explained = nan(numel(norm_labels), n_comps, n_cv);

methods = {'none', 'zscore', 'soft'};

rng(42);
for m = 1:numel(methods)

    % Normalize
    FR_norm = normalize_FR(FR, methods{m});

    % Collapse time into observations: reshape to [trials*bins x units]
    % Each time bin for each trial is one observation
    [nt, nu, nb] = size(FR_norm);
    X = reshape(permute(FR_norm, [1 3 2]), nt * nb, nu);

    % Trial indices for splitting (split at trial level to keep time intact)
    trial_perm = randperm(nt);
    half       = floor(nt / 2);

    for cv = 1:n_cv
        perm         = randperm(nt);
        idx_A        = perm(1:half);
        idx_B        = perm(half+1:end);

        % Get [obs x units] matrices for each half
        % obs = trials * bins
        rows_A = reshape((idx_A' - 1) * nb + (1:nb), 1, []);
        rows_B = reshape((idx_B' - 1) * nb + (1:nb), 1, []);

        X_A = X(rows_A, :);
        X_B = X(rows_B, :);

        % Mean-center each half separately
        mu_A = mean(X_A, 1);
        X_A  = X_A - mu_A;
        X_B  = X_B - mu_A;  % center B using A's mean

        % PCA on half A
        C_A      = (X_A' * X_A) / (size(X_A, 1) - 1);
        [V, ~]   = eig(C_A);
        [~, ord] = sort(diag(eig(C_A)), 'descend');
        V        = V(:, ord);
        V        = V(:, 1:n_comps);

        % Project both halves onto PCs from A
        proj_A = X_A * V;   % [obs_A x n_comps]
        proj_B = X_B * V;   % [obs_B x n_comps]

        % Cross-validated variance explained:
        % cov between projections of A and B onto each PC
        % Use the minimum of the two half sizes for pairing
        n_pair  = min(size(proj_A,1), size(proj_B,1));
        cv_cov  = nan(1, n_comps);
        for k = 1:n_comps
            a = proj_A(1:n_pair, k);
            b = proj_B(1:n_pair, k);
            c = cov(a, b);
            cv_cov(k) = c(1, 2);
        end

        % Negative values set to zero (no cross-validated variance)
        cv_cov = max(cv_cov, 0);

        % Normalize to get proportion of variance explained
        total = sum(cv_cov);
        if total > 0
            cv_var_explained(m, :, cv) = cv_cov / total;
        end
    end
end

% Average and cumulate across CV folds
mean_cv_var  = squeeze(nanmean(cv_var_explained, 3));   % [n_methods x n_comps]
cum_cv_var   = cumsum(mean_cv_var, 2);

% Effective dimensionality: participation ratio
% D = (sum lambda)^2 / sum(lambda^2)
eff_dim = nan(1, numel(methods));
for m = 1:numel(methods)
    lam        = mean_cv_var(m, :);
    eff_dim(m) = sum(lam)^2 / sum(lam.^2);
end

fprintf('\nEffective Dimensionality (Participation Ratio):\n');
for m = 1:numel(methods)
    fprintf('  %-20s: %.2f\n', norm_labels{m}, eff_dim(m));
end

% -------------------------------------------------------------------------
% Plot
% -------------------------------------------------------------------------
colors = [0.2  0.2  0.2;    % no norm  (dark)
          0.29 0.47 0.81;   % z-score  (blue)
          0.84 0.37 0.37];  % soft     (red)

n_show = min(20, n_comps);  % show first 20 PCs

figure;
set(gcf, 'Color', 'w', 'Position', [100 100 1200 450]);

% --- Panel 1: Scree plot (per-component CV variance explained) ---
subplot(1, 3, 1);
hold on;
for m = 1:numel(methods)
    plot(1:n_show, mean_cv_var(m, 1:n_show) * 100, ...
        'o-', 'Color', colors(m,:), 'LineWidth', 2, 'MarkerSize', 5, ...
        'DisplayName', norm_labels{m});
end
xlabel('Principal Component', 'FontSize', 11);
ylabel('CV Variance Explained (%)', 'FontSize', 11);
title('Scree Plot', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'northeast');
box off;

% --- Panel 2: Cumulative variance explained ---
subplot(1, 3, 2);
hold on;
for m = 1:numel(methods)
    plot(1:n_show, cum_cv_var(m, 1:n_show) * 100, ...
        'o-', 'Color', colors(m,:), 'LineWidth', 2, 'MarkerSize', 5, ...
        'DisplayName', norm_labels{m});
end
yline(80, 'k--', '80%', 'LineWidth', 1.2, 'LabelHorizontalAlignment', 'left');
yline(95, 'k:',  '95%', 'LineWidth', 1.2, 'LabelHorizontalAlignment', 'left');
xlabel('Number of Components', 'FontSize', 11);
ylabel('Cumulative CV Variance Explained (%)', 'FontSize', 11);
title('Cumulative Variance', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'southeast');
box off;

% --- Panel 3: Effective dimensionality bar chart ---
subplot(1, 3, 3);
hold on;
for m = 1:numel(methods)
    bar(m, eff_dim(m), 0.5, 'FaceColor', colors(m,:), ...
        'EdgeColor', 'none', 'FaceAlpha', 0.85);
    text(m, eff_dim(m) + 0.1, sprintf('%.2f', eff_dim(m)), ...
        'HorizontalAlignment', 'center', 'FontSize', 10);
end
xticks(1:numel(methods));
xticklabels(norm_labels);
ylabel('Effective Dimensionality', 'FontSize', 11);
title('Participation Ratio', 'FontSize', 12, 'FontWeight', 'bold');
xlim([0.5, numel(methods) + 0.5]);
box off;

sgtitle('Cross-Validated PCA: Dimensionality by Normalization Method', ...
        'FontSize', 13, 'FontWeight', 'bold');