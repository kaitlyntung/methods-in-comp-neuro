% We will do cross-validated PCA on M1 activity to determine the true low 
% dimensional structure of the neuron's responses, comparing dimensionality 
% under no normalization, z-scoring, and soft normalization.

function results = linearDimensionalityReduction(R, varargin)
p = inputParser;
addParameter(p, 'PreWindow', 150, @isnumeric);
addParameter(p, 'PostWindow', 50, @isnumeric);
addParameter(p, 'BinSize', 20, @isnumeric);
addParameter(p, 'nFolds', 5, @isnumeric);
addParameter(p, 'SoftNorm', 5, @isnumeric);
parse(p, varargin{:});
opt = p.Results;

pre_window  = opt.PreWindow;
post_window = opt.PostWindow;
bin_size = opt.BinSize;
k = opt.nFolds;
soft_norm = opt.SoftNorm;

conditionIDs = [R.conditionID];
straight_trials = find(mod(conditionIDs, 3) == 1);
curved_trials = find(mod(conditionIDs, 3) == 2);
all_trials = [straight_trials(:); curved_trials(:)];
n_trials = numel(all_trials);
n_units = numel(R(1).unit);
n_straight = numel(straight_trials);
n_curved = numel(curved_trials);

time_bins = -pre_window : bin_size : post_window;
bin_centers = time_bins(1:end-1) + bin_size/2;
n_bins = numel(bin_centers);

onset_align = [R.moveOnsetTime];
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
disp('Firing rates created.');

fr_straight = squeeze(mean(firing_rates(1:n_straight, :, :), 1)); 
fr_curved = squeeze(mean(firing_rates(n_straight+1:end, :, :), 1));
fr_cond_avg = permute(cat(3, fr_straight', fr_curved'), [3 1 2]); 
n_conds = 2;

cond_labels  = [ones(n_straight, 1); 2*ones(n_curved, 1)];
fold_ids_trials = zeros(n_trials, 1);
for cond = 1:2
    cond_mask = find(cond_labels == cond);
    n_cond = numel(cond_mask);
    perm = randperm(n_cond);
    cond_folds = mod(0:n_cond-1, k) + 1;
    fold_ids_trials(cond_mask) = cond_folds(perm);
end

methods = {'none', 'zscore', 'soft'};
normalization_labels = {'No Normalization', 'Z-Score', 'Soft Normalization'};
colors_method = [0 0 0; 0 0 1; 1 0 0];
n_methods = numel(methods);

n_X_fixed = floor(n_units / 2);
max_dims = n_X_fixed - 1;
n_plot_dims = min(n_X_fixed - 1, 2*n_bins - 1);

cv_loss = nan(n_methods, k, k, max_dims);
all_explained = nan(n_methods, n_plot_dims);
all_score = cell(n_methods, 1);
all_coeff = cell(n_methods, 1);

for m = 1:n_methods
    for fold_r = 1:k
        train_idx = find(fold_ids_trials ~= fold_r);
        test_idx = find(fold_ids_trials == fold_r);
        n_train = numel(train_idx);
        n_test = numel(test_idx);

        for fold_c = 1:k
            neuron_perm = randperm(n_units);
            x_neurons = neuron_perm(1:n_X_fixed);
            y_neurons = neuron_perm(n_X_fixed+1:end);
            n_Y = numel(y_neurons);

            X_train = reshape(permute(firing_rates(train_idx, x_neurons, :), [1 3 2]), n_train*n_bins, n_X_fixed);
            Y_train = reshape(permute(firing_rates(train_idx, y_neurons, :), [1 3 2]), n_train*n_bins, n_Y);
            X_test = reshape(permute(firing_rates(test_idx,  x_neurons, :), [1 3 2]), n_test*n_bins,  n_X_fixed);
            Y_test = reshape(permute(firing_rates(test_idx,  y_neurons, :), [1 3 2]), n_test*n_bins,  n_Y);

            [X_train_n, X_test_n] = ldr_normalize(X_train, X_test, methods{m}, soft_norm);

            mu_X = mean(X_train_n, 1);
            X_train_c = X_train_n - mu_X;
            X_test_c = X_test_n  - mu_X;
            mu_Y = mean(Y_train, 1);
            Y_train_c = Y_train - mu_Y;

            [coeff, score, ~] = pca(X_train_c, 'NumComponents', max_dims);

            for d = 1:max_dims
                S_train = score(:, 1:d);
                beta = S_train \ Y_train_c;
                S_test = X_test_c * coeff(:, 1:d);
                Y_hat = S_test * beta + mu_Y;
                cv_loss(m, fold_r, fold_c, d) = mean((Y_test - Y_hat).^2, 'all');
            end
        end
    end

    FR_all = reshape(permute(firing_rates, [1 3 2]), n_trials*n_bins, n_units);
    FR_all_tmp = FR_all;
    [FR_all_n, ~] = ldr_normalize(FR_all_tmp, FR_all_tmp, methods{m}, soft_norm);

    FR_all_n = reshape(FR_all_n, n_trials, n_bins, n_units);
    psth_s = squeeze(mean(FR_all_n(1:n_straight, :, :), 1));
    psth_c = squeeze(mean(FR_all_n(n_straight+1:end, :, :), 1));

    X_plot_final = [psth_s; psth_c];
    X_plot_final = X_plot_final - mean(X_plot_final, 1);

    [coeff_full, score_full, ~, ~, explained] = pca(X_plot_final);
    all_explained(m, :) = explained(1:n_plot_dims)';
    all_score{m} = score_full;
    all_coeff{m} = coeff_full;
end

mean_loss = squeeze(mean(cv_loss, [2 3]));

ldr_plot_cv_loss(mean_loss, max_dims, n_methods, normalization_labels, colors_method);
ldr_plot_scree(all_explained, n_plot_dims, n_methods, normalization_labels, colors_method);
ldr_plot_cumulative(all_explained, n_plot_dims, n_methods, normalization_labels, colors_method);

results.cv_loss = cv_loss;
results.mean_loss = mean_loss;
results.explained = all_explained;
results.score = all_score;
results.coeff = all_coeff;
results.bin_centers = bin_centers;
results.params = opt;
end 

function [X_n, X_test_n] = ldr_normalize(X_train, X_test, method, soft_norm)
    switch method
        case 'none'
            X_n = X_train;
            X_test_n = X_test;
        case 'zscore'
            mu = mean(X_train, 1);
            sig = std(X_train, 0, 1);
            sig(sig == 0) = 1;
            X_n = (X_train - mu) ./ sig;
            X_test_n = (X_test  - mu) ./ sig;
        case 'soft'
            rng = max(X_train, [], 1) - min(X_train, [], 1);
            X_n = X_train ./ (rng + soft_norm);
            X_test_n = X_test ./ (rng + soft_norm);
    end
end

function ldr_plot_cv_loss(mean_loss, max_dims, n_methods, norm_labels, colors)
    figure('Color', 'w');
    for m = 1:n_methods
        subplot(1, n_methods, m);
        plot(1:max_dims, mean_loss(m,:), '-o', ...
            'Color', colors(m,:), 'LineWidth', 2);
        xlabel('Number of PCs (k)');
        ylabel('MSE');
        title(norm_labels{m});
        grid on;
    end
    sgtitle('Cross-validated PCA loss by dimensionality');
end

function ldr_plot_scree(all_explained, n_plot_dims, n_methods, norm_labels, colors)
    n_plot = min(20, n_plot_dims);
    figure('Color', 'w'); hold on;
    for m = 1:n_methods
        plot(1:n_plot, all_explained(m, 1:n_plot), '-o', ...
            'Color', colors(m,:), 'LineWidth', 2);
    end
    xlabel('Principal component');
    ylabel('Variance explained (%)');
    title('Scree plot');
    legend(norm_labels, 'Location', 'northeast');
    grid on;
end

function ldr_plot_cumulative(all_explained, n_plot_dims, n_methods, norm_labels, colors)
    n_plot = min(20, n_plot_dims);
    figure('Color', 'w'); hold on;
    for m = 1:n_methods
        cum_exp = cumsum(all_explained(m, :));
        plot(1:n_plot, cum_exp(1:n_plot), '-o', ...
            'Color', colors(m,:), 'LineWidth', 2);
    end
    yline(80, 'k--', 'LineWidth', 1.5, 'Label', '80%');
    yline(95, 'k:',  'LineWidth', 1.5, 'Label', '95%');
    xlabel('Number of components');
    ylabel('Cumulative variance explained (%)');
    title('Cumulative variance');
    ylim([0 100]);
    legend(norm_labels, 'Location', 'southeast');
    grid on;
end