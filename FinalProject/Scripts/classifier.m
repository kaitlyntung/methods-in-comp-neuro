% -------------------------------------------------------------------------
% SVM and Logistic Regression classification of curved vs straight reaches
% from M1 neural activity during delay and movement periods separately
% With Monte Carlo shuffle test for significance
% -------------------------------------------------------------------------

% Time windows (ms relative to movement onset)
windows.delay    = [-500, 0];
windows.movement = [0,    500];
window_labels    = {'Delay Period', 'Movement Period'};
window_names     = {'delay', 'movement'};
n_windows        = 2;

onset_align = [R.moveOnsetTime];
n_units     = numel(R(1).unit);
bin_size    = 50;  % ms

% Labels: 0 = straight, 1 = curved
all_trials  = [straight_trials(:); curved_trials(:)];
labels      = [zeros(numel(straight_trials), 1); ...
               ones(numel(curved_trials),    1)];
n_trials    = numel(all_trials);

n_perm      = 1000;
n_cv_folds  = 5;    % k-fold cross-validation for classifier accuracy
rng(42);

% -------------------------------------------------------------------------
% Step 1: Build mean firing rate features per trial per window
% Feature matrix X: [n_trials x n_units] — mean FR in window
% -------------------------------------------------------------------------
function X = build_features(R, all_trials, onset_align, n_units, window, bin_size)
    n_trials = numel(all_trials);
    X        = nan(n_trials, n_units);
    time_bins = window(1) : bin_size : window(2);
    n_bins    = numel(time_bins) - 1;

    for t = 1:n_trials
        tr = all_trials(t);
        t0 = onset_align(tr);
        for u = 1:n_units
            spike_times = R(tr).unit(u).spikeTimes;
            rel_times   = spike_times - t0;
            count = 0;
            for b = 1:n_bins
                count = count + sum(rel_times >= time_bins(b) & ...
                                    rel_times <  time_bins(b+1));
            end
            X(t, u) = (count / n_bins) / (bin_size / 1000);  % mean FR in Hz
        end
    end
end

% -------------------------------------------------------------------------
% Step 2: Cross-validated classification accuracy
% -------------------------------------------------------------------------
function acc = classify_kfold(X, y, method, k)
    cv      = cvpartition(y, 'KFold', k, 'Stratify', true);
    correct = 0;
    for fold = 1:k
        X_train = X(cv.training(fold), :);
        y_train = y(cv.training(fold));
        X_test  = X(cv.test(fold),     :);
        y_test  = y(cv.test(fold));

        % Z-score using training set statistics
        mu      = mean(X_train, 1);
        sigma   = std(X_train,  0, 1);
        sigma(sigma == 0) = 1;
        X_train = (X_train - mu) ./ sigma;
        X_test  = (X_test  - mu) ./ sigma;

        switch method
            case 'svm'
                mdl  = fitcsvm(X_train, y_train, ...
                               'KernelFunction', 'linear', ...
                               'Standardize',    false, ...
                               'BoxConstraint',  1);
                pred = predict(mdl, X_test);

            case 'logreg'
                mdl  = fitclinear(X_train, y_train, ...
                                  'Learner',      'logistic', ...
                                  'Regularization','ridge', ...
                                  'Lambda',        1/size(X_train,1));
                pred = predict(mdl, X_test);
        end
        correct = correct + sum(pred == y_test);
    end
    acc = correct / numel(y);
end

% -------------------------------------------------------------------------
% Step 3: Run classifiers and Monte Carlo shuffle test
% -------------------------------------------------------------------------
methods       = {'svm', 'logreg'};
method_labels = {'SVM (Linear)', 'Logistic Regression'};
n_methods     = numel(methods);

% Results storage
obs_acc   = nan(n_methods, n_windows);   % observed accuracy
perm_acc  = nan(n_methods, n_windows, n_perm);
p_values  = nan(n_methods, n_windows);

for w = 1:n_windows
    fprintf('Window: %s\n', window_labels{w});

    % Build feature matrix for this window
    X = build_features(R, all_trials, onset_align, n_units, ...
                       [windows.(window_names{w})], bin_size);

    for m = 1:n_methods
        fprintf('  Method: %s\n', method_labels{m});

        % Observed accuracy
        obs_acc(m, w) = classify_kfold(X, labels, methods{m}, n_cv_folds);
        fprintf('    Observed accuracy: %.3f\n', obs_acc(m, w));

        % Monte Carlo shuffle
        parfor p = 1:n_perm
            shuf_labels      = labels(randperm(n_trials));
            perm_acc(m,w, p) = classify_kfold(X, shuf_labels, methods{m}, n_cv_folds);
        end

        % P-value: proportion of permutations >= observed
        p_values(m, w) = mean(squeeze(perm_acc(m, w, :)) >= obs_acc(m, w));
        fprintf('    p-value: %.4f\n', p_values(m, w));
    end
end

% -------------------------------------------------------------------------
% Print summary table
% -------------------------------------------------------------------------
fprintf('\n%s\n', repmat('-', 1, 65));
fprintf('%-22s %-20s %-10s %-10s\n', 'Method', 'Window', 'Accuracy', 'p-value');
fprintf('%s\n', repmat('-', 1, 65));
for w = 1:n_windows
    for m = 1:n_methods
        fprintf('%-22s %-20s %-10.3f %-10.4f\n', ...
            method_labels{m}, window_labels{w}, ...
            obs_acc(m,w), p_values(m,w));
    end
end
fprintf('%s\n', repmat('-', 1, 65));

% -------------------------------------------------------------------------
% Plot
% -------------------------------------------------------------------------
colors = [0.29 0.47 0.81;   % SVM
          0.84 0.37 0.37];  % LogReg

figure;
set(gcf, 'Color', 'w', 'Position', [100 100 1200 500]);

for w = 1:n_windows
    subplot(1, n_windows, w);
    hold on;

    null_data = squeeze(perm_acc(:, w, :));   % [n_methods x n_perm]

    % Plot null distributions
    for m = 1:n_methods
        histogram(null_data(m, :), 40, ...
            'Normalization', 'probability', ...
            'FaceColor',     colors(m, :), ...
            'EdgeColor',     'none', ...
            'FaceAlpha',     0.4, ...
            'DisplayName',   sprintf('%s null', method_labels{m}));
    end

    % Observed accuracy lines
    for m = 1:n_methods
        xline(obs_acc(m, w), '-', ...
            'Color',       colors(m, :), ...
            'LineWidth',   2.5, ...
            'DisplayName', sprintf('%s observed (p=%.3f)', ...
                           method_labels{m}, p_values(m, w)));
    end

    % Chance line
    xline(0.5, 'k--', 'Chance', 'LineWidth', 1.5, ...
          'LabelHorizontalAlignment', 'left');

    xlabel('Classification Accuracy', 'FontSize', 11);
    ylabel('Proportion', 'FontSize', 11);
    title(sprintf('Null Distribution\n%s', window_labels{w}), ...
          'FontSize', 12, 'FontWeight', 'bold');
    legend('Location', 'northwest', 'FontSize', 8);
    box off;
end

% Summary bar plot
figure;
set(gcf, 'Color', 'w', 'Position', [100 100 600 450]);
hold on;

bar_x      = [1 2 4 5];   % group delay and movement with a gap
bar_vals   = [obs_acc(1,1), obs_acc(2,1), obs_acc(1,2), obs_acc(2,2)];
bar_colors = [colors(1,:); colors(2,:); colors(1,:); colors(2,:)];

for i = 1:4
    bar(bar_x(i), bar_vals(i), 0.6, ...
        'FaceColor', bar_colors(i,:), 'EdgeColor', 'none', 'FaceAlpha', 0.85);
end

% P-value annotations
p_vals_ordered = [p_values(1,1), p_values(2,1), p_values(1,2), p_values(2,2)];
y_top = max(bar_vals) * 1.08;
for i = 1:4
    p = p_vals_ordered(i);
    if p < 0.001
        sig_str = '***';
    elseif p < 0.01
        sig_str = '**';
    elseif p < 0.05
        sig_str = '*';
    else
        sig_str = 'n.s.';
    end
    text(bar_x(i), bar_vals(i) + 0.01, sig_str, ...
        'HorizontalAlignment', 'center', 'FontSize', 12);
end

yline(0.5, 'k--', 'Chance', 'LineWidth', 1.5, ...
      'LabelHorizontalAlignment', 'left');

xticks(bar_x);
xticklabels({'SVM', 'LogReg', 'SVM', 'LogReg'});
ylabel('Cross-Validated Accuracy', 'FontSize', 11);
ylim([0.4, y_top + 0.05]);

% Window labels above groups
text(1.5, y_top + 0.03, 'Delay Period',    'HorizontalAlignment', 'center', ...
     'FontSize', 11, 'FontWeight', 'bold');
text(4.5, y_top + 0.03, 'Movement Period', 'HorizontalAlignment', 'center', ...
     'FontSize', 11, 'FontWeight', 'bold');

% Legend patches
patch(nan, nan, colors(1,:), 'DisplayName', 'SVM');
patch(nan, nan, colors(2,:), 'DisplayName', 'Logistic Regression');
legend('Location', 'southeast');
box off;

title('Classifier Accuracy: Curved vs Straight Reaches', ...
      'FontSize', 12, 'FontWeight', 'bold');