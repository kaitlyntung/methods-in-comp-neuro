% We will do cross-validated PCA on M1 activity to determine the true low 
% dimensional structure of the neuron's responses, comparing dimensionality 
% under no normalization, z-scoring, and soft normalization.

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
disp('firing rates created')


n_straight = numel(straight_trials);
n_curved   = numel(curved_trials);

% Average firing rates across trials per condition [units x bins]
psth_straight = squeeze(mean(firing_rates(1:n_straight, :, :), 1));
psth_curved   = squeeze(mean(firing_rates(n_straight+1:end, :, :), 1));

% Stack conditions along time: [units x (2*bins)]
X_psth = [psth_straight, psth_curved];

normalization_labels = {'No Normalization', 'Z-Score', 'Soft Normalization'};
methods = {'none', 'zscore', 'soft'};
soft_norm_alpha = 5;
n_plot = 10;
colors_method = [0 0 0; 0 0 0.8; 0.8 0 0];  % black, blue, red

% Store results for each method
all_explained = nan(numel(methods), size(X_psth, 2));
all_score     = cell(numel(methods), 1);

for m = 1:numel(methods)
    switch methods{m}
        case 'none'
            X_norm = X_psth;
        case 'zscore'
            mu    = mean(X_psth, 2);
            sigma = std(X_psth, 0, 2);
            sigma(sigma == 0) = 1;
            X_norm = (X_psth - mu) ./ sigma;
        case 'soft'
            fr_range = max(X_psth, [], 2) - min(X_psth, [], 2);
            X_norm   = X_psth ./ (fr_range + soft_norm_alpha);
    end

    % Mean center across time
    X_norm = X_norm - mean(X_norm, 2);

    % PCA
    [~, score, ~, ~, explained] = pca(X_norm');
    all_explained(m, 1:numel(explained)) = explained;
    all_score{m} = score;
end

% Scree Plot
figure;
hold on;
for m = 1:numel(methods)
    plot(1:n_plot, all_explained(m, 1:n_plot), '-o', ...
        'Color', colors_method(m,:), 'LineWidth', 2);
end
xlabel('Principal Component');
ylabel('Variance Explained (%)');
title('Scree Plot');
legend(normalization_labels, 'Location', 'northeast');

% Cumulative Variance
figure;
hold on;
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

figure;
pc_colors = nebula(3);
for m = 1:numel(methods)
    subplot(1, 3, m);
    hold on;
    score = all_score{m};
    for k = 1:3
        pc_straight = score(1:n_bins, k);
        pc_curved = score(n_bins+1:end, k);
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

%%
figure;
condition_names = {'Straight', 'Curved'};
cond_mod = [1, 2];

for cond = 1:2
    subplot(1, 2, cond);
    hold on;
    
    for d = 1:36
        target_condID = straight_condIDs(d) + (cond_mod(cond) - 1);
        trial_idx = find(conditionIDs == target_condID);
        
        if numel(trial_idx) < 3
            continue;
        end
        
        % Average across trials → one PSTH per direction×condition
        psth = squeeze(mean(firing_rates(trial_idx, :, :), 1));  % units x bins
        psth_norm = psth ./ (fr_range + 5);
        psth_norm = psth_norm - mean(psth_norm, 2);
        
        % Project onto shared PC space → one trajectory
        score = (coeff(:,1:3)' * psth_norm)';  % bins x 3
        
        % Smooth trajectory a little to reduce jagginess
        score = smoothdata(score, 1, 'gaussian', 5);
        
        % Plot trajectory
        plot3(score(:,1), score(:,2), score(:,3), '-', ...
            'Color', dir_colors(d,:), 'LineWidth', 2);
        
        % Start dot
        plot3(score(1,1), score(1,2), score(1,3), ...
            'o', 'Color', dir_colors(d,:), ...
            'MarkerFaceColor', dir_colors(d,:), 'MarkerSize', 5);
        
        % Movement onset diamond
        plot3(score(onset_bin,1), score(onset_bin,2), score(onset_bin,3), ...
            'd', 'Color', dir_colors(d,:), ...
            'MarkerFaceColor', dir_colors(d,:), 'MarkerSize', 7);
        
        % End square
        plot3(score(end,1), score(end,2), score(end,3), ...
            's', 'Color', dir_colors(d,:), ...
            'MarkerFaceColor', dir_colors(d,:), 'MarkerSize', 5);
    end
    
    xlabel(sprintf('PC1 (%.1f%%)', explained(1)));
    ylabel(sprintf('PC2 (%.1f%%)', explained(2)));
    zlabel(sprintf('PC3 (%.1f%%)', explained(3)));
    title(condition_names{cond});
    grid on;
    view(45, 25);
    axis tight;
end

% Match axis limits across both panels
subplot(1,2,1); ax1 = gca;
subplot(1,2,2); ax2 = gca;
all_lims = [ax1.XLim; ax2.XLim; ax1.YLim; ax2.YLim; ax1.ZLim; ax2.ZLim];
x_lim = [min(all_lims([1,2],1)), max(all_lims([1,2],2))];
y_lim = [min(all_lims([3,4],1)), max(all_lims([3,4],2))];
z_lim = [min(all_lims([5,6],1)), max(all_lims([5,6],2))];

%%
% Build angle-to-conditionID mapping
targetXY = reshape([R.targetXY], 2, [])';
angles = atan2d(targetXY(:,2), targetXY(:,1));
angleGroups = round(angles / 45) * 45;
angleGroups(angleGroups == 180) = -180;
uniqueGroups = unique(angleGroups);
nGroups = length(uniqueGroups);  % should be 7

% Color by angle group — use circular colormap so opposite dirs are distinct
dir_colors = hsv(nGroups);

% Angle labels for legend
angle_labels = arrayfun(@(a) sprintf('%d°', a), uniqueGroups, 'UniformOutput', false);

%%
targetXY = reshape([R.targetXY], 2, [])';angles = atan2d(targetXY(:,2), targetXY(:,1));angleGroups = round(angles / 45) * 45;angleGroups(angleGroups == 180) = -180;uniqueGroups = unique(angleGroups);nGroups = length(uniqueGroups);

figure;
condition_names = {'Straight', 'Curved'};
cond_mod = [1, 2];

for cond = 1:2
    subplot(1, 2, cond);
    hold on;

    for d = 1:nGroups
        % Find trials matching this angle group AND this condition type
        angle_mask = angleGroups == uniqueGroups(d);
        cond_mask  = mod(conditionIDs, 3) == cond_mod(cond);
        trial_idx  = find(angle_mask & cond_mask);

        if numel(trial_idx) < 3
            continue;
        end

        % Average across trials → one PSTH per direction×condition
        psth = squeeze(mean(firing_rates(trial_idx, :, :), 1));  % units x bins
        psth_norm = psth ./ (fr_range + 5);
        psth_norm = psth_norm - mean(psth_norm, 2);

        % Project onto shared PC space → one trajectory
        score = (coeff(:,1:3)' * psth_norm)';  % bins x 3

        % Smooth trajectory
        score = smoothdata(score, 1, 'gaussian', 5);

        % Plot trajectory
        plot3(score(:,1), score(:,2), score(:,3), '-', ...
            'Color', dir_colors(d,:), 'LineWidth', 2);

        % Start dot
        plot3(score(1,1), score(1,2), score(1,3), ...
            'o', 'Color', dir_colors(d,:), ...
            'MarkerFaceColor', dir_colors(d,:), 'MarkerSize', 5);

        % Movement onset diamond
        plot3(score(onset_bin,1), score(onset_bin,2), score(onset_bin,3), ...
            'd', 'Color', dir_colors(d,:), ...
            'MarkerFaceColor', dir_colors(d,:), 'MarkerSize', 7);

        % End square
        plot3(score(end,1), score(end,2), score(end,3), ...
            's', 'Color', dir_colors(d,:), ...
            'MarkerFaceColor', dir_colors(d,:), 'MarkerSize', 5);
    end

    xlabel(sprintf('PC1 (%.1f%%)', explained(1)));
    ylabel(sprintf('PC2 (%.1f%%)', explained(2)));
    zlabel(sprintf('PC3 (%.1f%%)', explained(3)));
    title(condition_names{cond});
    grid on;
    view(45, 25);
    axis tight;

    % Add direction legend to each panel
    legend(arrayfun(@(a) sprintf('%d°', a), uniqueGroups', 'UniformOutput', false), ...
        'Location', 'bestoutside', 'FontSize', 8);
end

% Match axis limits across both panels
subplot(1,2,1); ax1 = gca;
subplot(1,2,2); ax2 = gca;
all_lims = [ax1.XLim; ax2.XLim; ax1.YLim; ax2.YLim; ax1.ZLim; ax2.ZLim];
x_lim = [min(all_lims([1,2],1)), max(all_lims([1,2],2))];
y_lim = [min(all_lims([3,4],1)), max(all_lims([3,4],2))];
z_lim = [min(all_lims([5,6],1)), max(all_lims([5,6],2))];
set(ax1, 'XLim', x_lim, 'YLim', y_lim, 'ZLim', z_lim);
set(ax2, 'XLim', x_lim, 'YLim', y_lim, 'ZLim', z_lim);

sgtitle('Neural Trajectories — 7 directions × straight/curved');