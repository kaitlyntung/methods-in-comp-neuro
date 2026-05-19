conditionIDs = [R.conditionID];
straight_idx = mod(conditionIDs, 3) == 1;
curved_idx = mod(conditionIDs, 3) == 2;
straight_trials = find(straight_idx);
curved_trials = find(curved_idx);

pre_window = 500;
post_window = 1500;
bin_size = 20;
time_bins = -pre_window : bin_size : post_window;
bin_centers = time_bins(1:end-1) + bin_size/2;
n_bins = numel(bin_centers);
n_units = numel(R(1).unit);

onset_align = [R.moveOnsetTime];
all_trials = [straight_trials(:); curved_trials(:)];
n_trials = numel(all_trials);

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

fr_straight = squeeze(mean(firing_rates(1:numel(straight_trials), :, :), 1));  % [n_units x n_bins]
fr_curved   = squeeze(mean(firing_rates(numel(straight_trials)+1:end, :, :), 1));
fr_cond_avg = permute(cat(3, fr_straight', fr_curved'), [3 1 2]);  % [2 x n_bins x n_units]
n_conds = 2;

k = 5;
soft_norm = 5;
methods = {'none', 'zscore', 'soft'};
normalization_labels = {'No Normalization', 'Z-Score', 'Soft Normalization'};
colors_method = [0 0 0; 0 0 1; 1 0 0];

n_straight = numel(straight_trials);
n_curved = numel(curved_trials);
cond_labels = [ones(n_straight, 1); 2*ones(n_curved, 1)];

fold_ids_trials = zeros(n_trials, 1);
for cond = 1:2
    cond_mask = find(cond_labels == cond);
    n_cond = numel(cond_mask);
    perm = randperm(n_cond);
    cond_folds = mod(0:n_cond-1, k) + 1;
    fold_ids_trials(cond_mask) = cond_folds(perm);
end

n_X_fixed = floor(n_units / 2);
max_dims = n_X_fixed - 1;

cv_loss = nan(numel(methods), k, k, max_dims);
n_plot_dims = min(n_X_fixed - 1, 2*n_bins - 1);
all_explained = nan(numel(methods), n_plot_dims);
all_score = cell(numel(methods), 1);
all_coeff = cell(numel(methods), 1);

for m = 1:numel(methods)

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

            switch methods{m}
                case 'none'
                    X_train_n = X_train; X_test_n = X_test;
                case 'zscore'
                    mu_x = mean(X_train, 1);
                    sig_x = std(X_train, 0, 1); sig_x(sig_x==0) = 1;
                    X_train_n = (X_train - mu_x) ./ sig_x;
                    X_test_n  = (X_test  - mu_x) ./ sig_x;
                case 'soft'
                    rng_x = max(X_train,[],1) - min(X_train,[],1);
                    X_train_n = X_train ./ (rng_x + soft_norm);
                    X_test_n  = X_test  ./ (rng_x + soft_norm);
            end

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

    X_plot = reshape(fr_cond_avg(:, :, :), n_conds*n_bins, n_units);

    % switch methods{m}
    %     case 'none'
    %         X_plot_n = X_plot;
    %     case 'zscore'
    %         mu_x = mean(X_plot, 1);
    %         sig_x = std(X_plot, 0, 1); sig_x(sig_x==0) = 1;
    %         X_plot_n = (X_plot - mu_x) ./ sig_x;
    %     case 'soft'
    %         rng_x = max(X_plot,[],1) - min(X_plot,[],1);
    %         X_plot_n = X_plot ./ (rng_x + soft_norm);
    % end
    % Reshape single-trial firing rates: [n_trials*n_bins x n_units]
    FR_all = reshape(permute(firing_rates, [1 3 2]), n_trials*n_bins, n_units);
    
    % Normalize using same method as CV, but fit on all data for plotting
    switch methods{m}
        case 'none'
            FR_all_n = FR_all;
        case 'zscore'
            mu_x = mean(FR_all, 1);
            sig_x = std(FR_all, 0, 1); sig_x(sig_x==0) = 1;
            FR_all_n = (FR_all - mu_x) ./ sig_x;
        case 'soft'
            rng_x = max(FR_all,[],1) - min(FR_all,[],1);
            FR_all_n = FR_all ./ (rng_x + soft_norm);
    end
    
    % Reshape back to [n_trials x n_bins x n_units], average per condition
    FR_all_n = reshape(FR_all_n, n_trials, n_bins, n_units);
    psth_s = squeeze(mean(FR_all_n(1:n_straight, :, :), 1));      % [n_bins x n_units]
    psth_c = squeeze(mean(FR_all_n(n_straight+1:end, :, :), 1));  % [n_bins x n_units]
    
    % Stack cleanly: straight first, then curved
    X_plot_final = [psth_s; psth_c];                   % [2*n_bins x n_units]
    X_plot_final = X_plot_final - mean(X_plot_final, 1);
    
    [coeff_full, score_full, ~, ~, explained] = pca(X_plot_final);
    all_score{m} = score_full;  % rows 1:n_bins = straight, n_bins+1:end = curved
    % 
    % X_plot_n = X_plot_n - mean(X_plot_n, 1);
    % [coeff_full, score_full, ~, ~, explained] = pca(X_plot_n);

    all_explained(m, :) = explained(1:n_plot_dims)';
    all_score{m} = score_full;
    all_coeff{m} = coeff_full;
end

mean_loss = squeeze(mean(cv_loss, [2 3]));

%% CV loss plot
figure;
for m = 1:numel(methods)
    subplot(1, 3, m);
    plot(1:max_dims, mean_loss(m, :), '-o', 'Color', colors_method(m,:), 'LineWidth', 2);
    xlabel('Number of PCs (k)');
    ylabel('MSE');
    title(normalization_labels{m});
    grid on;
end
sgtitle('Cross-validated PCA loss by dimensionality');

%% Scree plot
n_plot = min(20, n_plot_dims);
figure; hold on;
for m = 1:numel(methods)
    plot(1:n_plot, all_explained(m, 1:n_plot), '-o', 'Color', colors_method(m,:), 'LineWidth', 2);
end
xlabel('Principal component');
ylabel('Variance explained (%)');
title('Scree plot');
legend(normalization_labels, 'Location', 'northeast');
grid on;

%% Cumulative variance
figure; hold on;
for m = 1:numel(methods)
    cum_exp = cumsum(all_explained(m, :));
    plot(1:n_plot, cum_exp(1:n_plot), '-o', 'Color', colors_method(m,:), 'LineWidth', 2);
end
yline(80, 'k--', 'LineWidth', 1.5, 'Label', '80%');
yline(95, 'k:',  'LineWidth', 1.5, 'Label', '95%');
xlabel('Number of components');
ylabel('Cumulative variance explained (%)');
title('Cumulative variance');
ylim([0 100]);
legend(normalization_labels, 'Location', 'southeast');
grid on;

%% 2D PC trajectories
figure;
pc_colors = nebula(3);
for m = 1:numel(methods)
    subplot(1, 3, m); hold on;
    score = all_score{m};
    pc_straight = score(1:n_bins, pc_idx);
    pc_curved   = score(n_bins+1:end, pc_idx);
    for pc = 1:3
        plot(bin_centers, score(1:n_bins, pc), '-',  'Color', pc_colors(pc,:), 'LineWidth', 2);
        plot(bin_centers, score(n_bins+1:end, pc), '--', 'Color', pc_colors(pc,:), 'LineWidth', 2);
    end
    xline(0, 'k--', 'LineWidth', 1.5);
    xlabel('Time from movement onset (ms)');
    ylabel('PC score');
    title(normalization_labels{m});
    if m == 1
        legend('PC1 straight','PC1 curved','PC2 straight','PC2 curved', ...
               'PC3 straight','PC3 curved','Location','best');
    end
end
sgtitle('Top 3 PCs Over Time');
