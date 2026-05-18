% We will classify curved vs. straight reaches based on neural data from
% the delay and movement period separately using SVM and logistic
% regression to assess when kinematic information is encoded in M1.
% We will assess the significance of our classifier using a Monte Carlo
% shuffle test, where we shuffle curved vs straight trial labels 1000 times
% to generate a null distribution and compute a p-value for each
% method/time-window combination.

windows.delay    = [-500, 0];
windows.movement = [0, 500];
window_labels    = {'Delay Period', 'Movement Period'};
window_names     = {'delay', 'movement'};
n_windows        = numel(window_labels);   % <-- was missing

all_trials = [straight_trials(:); curved_trials(:)];
labels     = [zeros(numel(straight_trials), 1); ones(numel(curved_trials), 1)];
n_trials   = numel(all_trials);

n_perm     = 1000;
n_cv_folds = 3;

delay_bins    = bin_centers >= -500 & bin_centers < 0;
movement_bins = bin_centers >= 0    & bin_centers < 500;

% firing rates from the linearDimensionalityReduction script
X_delay    = squeeze(mean(firing_rates(:, :, delay_bins),    3));
X_movement = squeeze(mean(firing_rates(:, :, movement_bins), 3));

window_features = {X_delay, X_movement};

methods       = {'svm', 'logreg'};
method_labels = {'SVM (Linear)', 'Logistic Regression'};
n_methods     = numel(methods);

obs_acc  = nan(n_methods, n_windows);
perm_acc = nan(n_methods, n_windows, n_perm);
p_values = nan(n_methods, n_windows);

for w = 1:n_windows
    X = window_features{w};
    for m = 1:n_methods
        obs_acc(m, w) = classify_kfold(X, labels, methods{m}, n_cv_folds);

        % parfor moved to the outermost loop it can own cleanly;
        % broadcast X, labels, method string, n_cv_folds as sliced/const inputs
        method_m = methods{m};
        perm_acc_tmp = nan(n_perm, 1);

        parfor p = 1:n_perm
            disp(p)
            shuf_labels      = labels(randperm(n_trials));
            perm_acc_tmp(p)  = classify_kfold(X, shuf_labels, method_m, n_cv_folds);
        end

        perm_acc(m, w, :) = perm_acc_tmp;
        p_values(m, w)    = (sum(perm_acc_tmp >= obs_acc(m, w)) + 1) / (n_perm + 1);

        fprintf('Done: %s | %s\n', method_labels{m}, window_labels{w});
    end
end

% --- Results table ---
fprintf('\n%s\n', repmat('-', 1, 65));
fprintf('%-22s %-20s %-10s %-10s\n', 'Method', 'Window', 'Accuracy', 'p-value');
fprintf('%s\n', repmat('-', 1, 65));
for w = 1:n_windows
    for m = 1:n_methods
        fprintf('%-22s %-20s %-10.3f %-10.4f %s\n', ...
            method_labels{m}, window_labels{w}, ...
            obs_acc(m, w), p_values(m, w), significance_label(p_values(m, w)));
    end
end
fprintf('%s\n', repmat('-', 1, 65));


%% PLOT
clr.svm    = [0.20, 0.45, 0.75];
clr.logreg = [0.85, 0.33, 0.10];
method_colors = {clr.svm, clr.logreg};

chance_level = 0.50;

% --- Figure 1: Accuracy comparison ---
figure('Name', 'Accuracy Comparison', 'Color', 'w', 'Position', [100, 100, 560, 420]);

group_centers = 1:n_windows;
bar_width     = 0.3;
offsets       = [-0.18, 0.18];

ax = axes; hold(ax, 'on');
bar_handles = gobjects(n_methods, 1);

for m = 1:n_methods
    x_pos = group_centers + offsets(m);
    bar_handles(m) = bar(x_pos, obs_acc(m, :), bar_width, ...
                         'FaceColor', method_colors{m}, ...
                         'EdgeColor', 'none', ...
                         'FaceAlpha', 0.85);

    null_95 = squeeze(prctile(perm_acc(m, :, :), 95, 3));
    for w = 1:n_windows
        plot(ax, [x_pos(w), x_pos(w)], [obs_acc(m,w), null_95(w)], ...
             'k-', 'LineWidth', 1.2);
        plot(ax, x_pos(w) + [-0.05, 0.05], [null_95(w), null_95(w)], ...
             'k-', 'LineWidth', 1.2);

        star  = significance_label(p_values(m, w));
        y_top = max(obs_acc(m,w), null_95(w)) + 0.025;

        text(ax, x_pos(w), obs_acc(m,w) - 0.03, ...
             sprintf('%.1f%%', obs_acc(m,w) * 100), ...
             'HorizontalAlignment', 'center', ...
             'FontSize', 9, 'Color', 'w', 'FontWeight', 'bold');

        if ~strcmp(star, '(n.s.)')
            text(ax, x_pos(w), y_top, star, ...
                 'HorizontalAlignment', 'center', ...
                 'FontSize', 13, 'FontWeight', 'bold', 'Color', 'k');
        end
    end
end

yline(ax, chance_level, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.4, ...
      'Label', 'Chance (50%)', 'LabelHorizontalAlignment', 'left', 'FontSize', 10);

ax.XTick         = group_centers;
ax.XTickLabel    = window_labels;
ax.XLim          = [0.5, n_windows + 0.5];
ax.YLim          = [0.3, 1.05];
ax.YLabel.String = 'Accuracy';
ax.Title.String  = 'Classifier Accuracy by Method & Time Window';
ax.FontSize      = 12;
ax.Box           = 'off';
legend(ax, bar_handles, method_labels, 'Location', 'northwest', 'Box', 'off');

% --- Figure 2: Null distributions ---
figure('Name', 'Null Distributions', 'Color', 'w', ...
       'Position', [150, 150, 320 * n_windows, 280 * n_methods]);

panel = 0;
for m = 1:n_methods
    for w = 1:n_windows
        panel = panel + 1;
        ax2   = subplot(n_methods, n_windows, panel);
        hold(ax2, 'on');

        null_dist = squeeze(perm_acc(m, w, :));

        histogram(null_dist, 30, ...
                  'FaceColor', method_colors{m}, 'EdgeColor', 'none', ...
                  'FaceAlpha', 0.55, 'Normalization', 'probability');

        thresh = prctile(null_dist, 95);
        xline(ax2, thresh, '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.4, ...
              'Label', '95th pct.', 'LabelVerticalAlignment', 'bottom', 'FontSize', 8);

        xline(ax2, obs_acc(m, w), '-', 'Color', method_colors{m}, 'LineWidth', 2.2, ...
              'Label', sprintf('Observed (%.1f%%)', obs_acc(m,w) * 100), ...
              'LabelVerticalAlignment', 'top', 'FontSize', 8);

        xline(ax2, chance_level, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2, ...
              'Label', 'Chance', 'LabelVerticalAlignment', 'bottom', 'FontSize', 8);

        % --- p-value formatted to never display as 0.0000 ---
        p_val = p_values(m, w);
        if p_val < 0.0001
            p_str = '< 0.0001';
        else
            p_str = sprintf('%.4f', p_val);
        end

        ax2.XLabel.String = 'Accuracy';
        ax2.YLabel.String = 'Proportion';
        ax2.Title.String  = sprintf('%s — %s\nacc = %.1f%%,  p = %s  %s', ...
            method_labels{m}, window_labels{w}, ...
            obs_acc(m,w) * 100, p_str, significance_label(p_val));
        ax2.FontSize = 10;
        ax2.Box      = 'off';
    end
end
sgtitle('Null Distributions vs. Observed Accuracy', 'FontSize', 13, 'FontWeight', 'bold');


%% Helper functions
function acc = classify_kfold(X, y, method, k)
    cv      = cvpartition(y, 'KFold', k, 'Stratify', true);
    correct = 0;
    for fold = 1:k
        X_train = X(cv.training(fold), :);
        y_train = y(cv.training(fold));
        X_test  = X(cv.test(fold),     :);
        y_test  = y(cv.test(fold));

        mu             = mean(X_train, 1);
        sigma          = std(X_train,  0, 1);
        sigma(sigma==0)= 1;
        X_train        = (X_train - mu) ./ sigma;
        X_test         = (X_test  - mu) ./ sigma;

        switch method
            case 'svm'
                mdl = fitcsvm(X_train, y_train, ...
                              'KernelFunction', 'linear', ...
                              'Standardize',    false, ...
                              'BoxConstraint',  1);
            case 'logreg'
                mdl = fitclinear(X_train, y_train, ...
                                 'Learner',        'logistic', ...
                                 'Regularization', 'ridge', ...
                                 'Lambda',         1/size(X_train,1));
        end
        pred    = predict(mdl, X_test);
        correct = correct + sum(pred == y_test);
    end
    acc = correct / numel(y);
end

function label = significance_label(p)
    if     p < 0.001;  label = '***';
    elseif p < 0.01;   label = '**';
    elseif p < 0.05;   label = '*';
    else;              label = '(n.s.)';
    end
end