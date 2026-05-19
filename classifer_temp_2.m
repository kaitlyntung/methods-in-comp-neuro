
%% =========================================================================
% DEFINE STRAIGHT VS CURVED TRIALS
% =========================================================================

conditionIDs = [R.conditionID];

straight_idx = mod(conditionIDs, 3) == 1;
curved_idx   = mod(conditionIDs, 3) == 2;

straight_trials = find(straight_idx);
curved_trials   = find(curved_idx);

%% =========================================================================
% DEFINE TIME WINDOWS
% =========================================================================

windows.delay    = [-500, 0];
windows.movement = [0, 500];

window_labels = {'Delay Period', 'Movement Period'};
window_names  = {'delay', 'movement'};

n_windows = numel(window_labels);

%% =========================================================================
% DEFINE TIME BINS
% =========================================================================

bin_size = 50;

time_bins = windows.delay(1):bin_size:windows.movement(2);

bin_centers = time_bins(1:end-1) + bin_size/2;

delay_bins = ...
    bin_centers >= -500 & ...
    bin_centers < 0;

movement_bins = ...
    bin_centers >= 0 & ...
    bin_centers < 500;

%% =========================================================================
% DEFINE CLASSIFIERS
% =========================================================================

methods = {'svm', 'logreg'};

method_labels = { ...
    'SVM (Linear)', ...
    'Logistic Regression'};

n_methods = numel(methods);

%% =========================================================================
% DEFINE REACH DIRECTION LABELS
% =========================================================================
%
% Convert target locations into direction categories.
%
% Example:
%   -180°, -135°, -90°, ...
%
% =========================================================================

targetXY = reshape([R.targetXY], 2, [])';

angles = atan2d(targetXY(:,2), targetXY(:,1));

angleGroups = round(angles / 45) * 45;

angleGroups(angleGroups == 180) = -180;

uniqueGroups = unique(angleGroups);

nGroups = numel(uniqueGroups);

% Integer labels for classifiers
dir_labels_all = nan(numel(R),1);

for g = 1:nGroups

    dir_labels_all(angleGroups == uniqueGroups(g)) = g;

end

%% =========================================================================
% CHANCE LEVEL
% =========================================================================

chance_level_dir = 1 / nGroups;

%% =========================================================================
% DEFINE CONDITIONS
% =========================================================================

barrier_conditions = { ...
    straight_trials, ...
    curved_trials};

barrier_labels_str = { ...
    'Straight (No Barriers)', ...
    'Curved (Barriers)'};

n_barrier = numel(barrier_conditions);

%% =========================================================================
% ANALYSIS SETTINGS
% =========================================================================

n_perm     = 1000;
n_cv_folds = 5;

%% =========================================================================
% PREALLOCATE STORAGE
% =========================================================================

obs_acc_dir  = nan(n_methods, n_windows, n_barrier);

perm_acc_dir = nan(n_methods, n_windows, n_barrier, n_perm);

p_values_dir = nan(n_methods, n_windows, n_barrier);

pred_store = cell(n_methods, n_windows, n_barrier);

true_store = cell(n_methods, n_windows, n_barrier);

%% =========================================================================
% MAIN ANALYSIS LOOP
% =========================================================================

for b = 1:n_barrier

    %% --------------------------------------------------------------------
    % Select trials for this condition
    % ---------------------------------------------------------------------

    trial_ids = barrier_conditions{b};

    dir_labels = dir_labels_all(trial_ids);

    n_tr = numel(trial_ids);

    %% --------------------------------------------------------------------
    % Extract firing rates for THIS subset of trials
    % ---------------------------------------------------------------------

    X_delay_b = squeeze(mean( ...
        firing_rates(trial_ids,:,delay_bins), 3));

    X_movement_b = squeeze(mean( ...
        firing_rates(trial_ids,:,movement_bins), 3));

    window_features_b = { ...
        X_delay_b, ...
        X_movement_b};

    %% --------------------------------------------------------------------
    % CREATE FIXED CROSS-VALIDATION PARTITION
    %
    % IMPORTANT:
    % SAME partition reused for:
    %   - observed decoding
    %   - all permutations
    %
    % ---------------------------------------------------------------------

    cv = cvpartition( ...
        dir_labels, ...
        'KFold', n_cv_folds, ...
        'Stratify', true);

    %% --------------------------------------------------------------------
    % Loop over time windows
    % ---------------------------------------------------------------------

    for w = 1:n_windows

        X = window_features_b{w};

        %% ----------------------------------------------------------------
        % Loop over classifier methods
        % -----------------------------------------------------------------

        for m = 1:n_methods

            fprintf('\n========================================\n');
            fprintf('%s | %s | %s\n', ...
                method_labels{m}, ...
                window_labels{w}, ...
                barrier_labels_str{b});
            fprintf('========================================\n');

            %% ------------------------------------------------------------
            % OBSERVED DECODING
            % -------------------------------------------------------------

            [pred, true_lab] = classify_kfold_predictions( ...
                X, ...
                dir_labels, ...
                methods{m}, ...
                cv);

            pred_store{m,w,b} = pred;
            true_store{m,w,b} = true_lab;

            obs_acc_dir(m,w,b) = mean(pred == true_lab);

            fprintf('Observed accuracy = %.3f\n', ...
                obs_acc_dir(m,w,b));

            %% ------------------------------------------------------------
            % MONTE CARLO PERMUTATION TEST
            % -------------------------------------------------------------
            %
            % Shuffle labels while keeping:
            %   - neural data fixed
            %   - CV folds fixed
            %
            % This estimates the NULL distribution:
            %
            % "What accuracies occur by chance?"
            %
            % -------------------------------------------------------------

            method_m = methods{m};

            perm_acc_tmp = nan(n_perm,1);

            parfor p = 1:n_perm

                % Shuffle labels
                shuf_labels = dir_labels(randperm(n_tr));

                % Decode using SAME CV partition
                perm_acc_tmp(p) = classify_kfold( ...
                    X, ...
                    shuf_labels, ...
                    method_m, ...
                    cv);

            end

            %% ------------------------------------------------------------
            % Store null distribution
            % -------------------------------------------------------------

            perm_acc_dir(m,w,b,:) = perm_acc_tmp;

            %% ------------------------------------------------------------
            % COMPUTE P-VALUE
            % -------------------------------------------------------------
            %
            % Monte Carlo permutation p-value:
            %
            % p =
            % (# permuted accuracies >= observed + 1)
            % ---------------------------------------
            %           (n_perm + 1)
            %
            % -------------------------------------------------------------

            p_values_dir(m,w,b) = ...
                (sum(perm_acc_tmp >= obs_acc_dir(m,w,b)) + 1) ...
                / (n_perm + 1);

            fprintf('p-value = %.4f\n', ...
                p_values_dir(m,w,b));

        end
    end
end

%% =========================================================================
% RESULTS TABLE
% =========================================================================

fprintf('\n%s\n', repmat('-',1,75));

fprintf('%-22s %-20s %-24s %-10s %-10s\n', ...
    'Method', ...
    'Window', ...
    'Condition', ...
    'Accuracy', ...
    'p-value');

fprintf('%s\n', repmat('-',1,75));

for b = 1:n_barrier

    for w = 1:n_windows

        for m = 1:n_methods

            p_val = p_values_dir(m,w,b);

            if p_val < 0.0001

                p_str = '< 0.0001';

            else

                p_str = sprintf('%.4f', p_val);

            end

            fprintf('%-22s %-20s %-24s %-10.3f %-10s %s\n', ...
                method_labels{m}, ...
                window_labels{w}, ...
                barrier_labels_str{b}, ...
                obs_acc_dir(m,w,b), ...
                p_str, ...
                significance_label(p_val));

        end
    end
end

fprintf('%s\n', repmat('-',1,75));

%% =========================================================================
% CONFUSION MATRIX PLOTS
% =========================================================================

angle_strs = arrayfun( ...
    @(a) sprintf('%d°', a), ...
    uniqueGroups, ...
    'UniformOutput', false);

for b = 1:n_barrier

    figure( ...
        'Color', 'w', ...
        'Position', [100,100,420*n_windows,380*n_methods], ...
        'Name', barrier_labels_str{b});

    panel = 0;

    for m = 1:n_methods

        for w = 1:n_windows

            panel = panel + 1;

            ax = subplot(n_methods, n_windows, panel);

            pred     = pred_store{m,w,b};
            true_lab = true_store{m,w,b};

            %% ------------------------------------------------------------
            % Compute confusion matrix
            % -------------------------------------------------------------

            C = zeros(nGroups, nGroups);

            for i = 1:nGroups

                for j = 1:nGroups

                    C(i,j) = sum( ...
                        true_lab == i & pred == j);

                end
            end

            %% ------------------------------------------------------------
            % Row-normalize
            % -------------------------------------------------------------

            C_norm = C ./ sum(C,2);

            %% ------------------------------------------------------------
            % Plot matrix
            % -------------------------------------------------------------

            imagesc(ax, C_norm);

            colormap(ax, 'sky');

            clim(ax, [0 1]);

            hold(ax, 'on');

            %% ------------------------------------------------------------
            % Add percentages inside cells
            % -------------------------------------------------------------

            for i = 1:nGroups

                for j = 1:nGroups

                    val = C_norm(i,j);

                    txt_col = 'k';

                    if val > 0.6
                        txt_col = 'w';
                    end

                    text(ax, j, i, ...
                        sprintf('%.0f%%', val*100), ...
                        'HorizontalAlignment', 'center', ...
                        'VerticalAlignment', 'middle', ...
                        'FontSize', 8, ...
                        'Color', txt_col);

                end
            end

            %% ------------------------------------------------------------
            % Highlight diagonal
            % -------------------------------------------------------------

            for i = 1:nGroups

                rectangle(ax, ...
                    'Position', [i-0.5, i-0.5, 1, 1], ...
                    'EdgeColor', [0.9 0.5 0.1], ...
                    'LineWidth', 1.5);

            end

            %% ------------------------------------------------------------
            % Title information
            % -------------------------------------------------------------

            p_val = p_values_dir(m,w,b);

            if p_val < 0.0001

                p_str = '< 0.0001';

            else

                p_str = sprintf('%.4f', p_val);

            end

            overall_acc = mean(pred == true_lab) * 100;

            ax.XTick = 1:nGroups;
            ax.YTick = 1:nGroups;

            ax.XTickLabel = angle_strs;
            ax.YTickLabel = angle_strs;

            ax.TickLength = [0 0];

            ax.FontSize = 9;

            xlabel(ax, 'Predicted direction');

            ylabel(ax, 'True direction');

            title(ax, sprintf( ...
                '%s | %s\nacc = %.1f%%, p = %s %s', ...
                method_labels{m}, ...
                window_labels{w}, ...
                overall_acc, ...
                p_str, ...
                significance_label(p_val)));

            cb = colorbar(ax);

            cb.Label.String = 'Proportion';

            cb.FontSize = 8;

        end
    end

    sgtitle(sprintf( ...
        'Direction Confusion Matrix — %s', ...
        barrier_labels_str{b}), ...
        'FontSize', 13, ...
        'FontWeight', 'bold');

end

%% =========================================================================
% HELPER FUNCTION:
% CROSS-VALIDATED DECODING ACCURACY
% =========================================================================

function acc = classify_kfold(X, y, method, cv)

    correct = 0;

    n_folds = cv.NumTestSets;

    for fold = 1:n_folds

        %% ------------------------------------------------------------
        % Train/test split
        % -------------------------------------------------------------

        train_idx = training(cv, fold);

        test_idx = test(cv, fold);

        X_train = X(train_idx,:);
        y_train = y(train_idx);

        X_test = X(test_idx,:);
        y_test = y(test_idx);

        %% ------------------------------------------------------------
        % Z-score normalization
        %
        % IMPORTANT:
        % Compute normalization ONLY on training data
        % -------------------------------------------------------------

        mu = mean(X_train,1);

        sigma = std(X_train,0,1);

        sigma(sigma == 0) = 1;

        X_train = (X_train - mu) ./ sigma;

        X_test = (X_test - mu) ./ sigma;

        %% ------------------------------------------------------------
        % Train classifier
        % -------------------------------------------------------------

        switch method

            case 'svm'

                t = templateSVM( ...
                    'KernelFunction', 'linear', ...
                    'Standardize', false, ...
                    'BoxConstraint', 1);

                mdl = fitcecoc( ...
                    X_train, ...
                    y_train, ...
                    'Learners', t);

            case 'logreg'

                t = templateLinear( ...
                    'Learner', 'logistic', ...
                    'Regularization', 'ridge', ...
                    'Lambda', 1/size(X_train,1));

                mdl = fitcecoc( ...
                    X_train, ...
                    y_train, ...
                    'Learners', t);

        end

        %% ------------------------------------------------------------
        % Predict test labels
        % -------------------------------------------------------------

        pred = predict(mdl, X_test);

        correct = correct + sum(pred == y_test);

    end

    %% ------------------------------------------------------------
    % Final decoding accuracy
    % -------------------------------------------------------------

    acc = correct / numel(y);

end

%% =========================================================================
% HELPER FUNCTION:
% STORE CROSS-VALIDATED PREDICTIONS
% =========================================================================

function [pred_all, true_all] = classify_kfold_predictions( ...
    X, y, method, cv)

    pred_all = nan(size(y));

    true_all = y;

    n_folds = cv.NumTestSets;

    for fold = 1:n_folds

        train_idx = training(cv, fold);

        test_idx = test(cv, fold);

        X_train = X(train_idx,:);
        y_train = y(train_idx);

        X_test = X(test_idx,:);

        %% ------------------------------------------------------------
        % Z-score normalization
        % -------------------------------------------------------------

        mu = mean(X_train,1);

        sigma = std(X_train,0,1);

        sigma(sigma == 0) = 1;

        X_train = (X_train - mu) ./ sigma;

        X_test = (X_test - mu) ./ sigma;

        %% ------------------------------------------------------------
        % Train classifier
        % -------------------------------------------------------------

        switch method

            case 'svm'

                t = templateSVM( ...
                    'KernelFunction', 'linear', ...
                    'Standardize', false, ...
                    'BoxConstraint', 1);

                mdl = fitcecoc( ...
                    X_train, ...
                    y_train, ...
                    'Learners', t);

            case 'logreg'

                t = templateLinear( ...
                    'Learner', 'logistic', ...
                    'Regularization', 'ridge', ...
                    'Lambda', 1/size(X_train,1));

                mdl = fitcecoc( ...
                    X_train, ...
                    y_train, ...
                    'Learners', t);

        end

        %% ------------------------------------------------------------
        % Store predictions
        % -------------------------------------------------------------

        pred_all(test_idx) = predict(mdl, X_test);

    end
end

%% =========================================================================
% HELPER FUNCTION:
% SIGNIFICANCE LABELS
% =========================================================================

function label = significance_label(p)

    if p < 0.001

        label = '***';

    elseif p < 0.01

        label = '**';

    elseif p < 0.05

        label = '*';

    else

        label = '(n.s.)';

    end
end