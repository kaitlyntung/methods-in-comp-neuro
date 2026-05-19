% We will classify curved vs. straight reaches based on neural data from 
% the delay and movement period separately using SVM and logistic 
% regression to assess when kinematic information is encoded in M1.

% We will assess the significance of our classifier using a Monte Carlo 
% shuffle test, where we shuffle curved vs straight trial labels 1000 times
% to generate a null distribution and compute a p-value for each
% method/time-window combination.

conditionIDs = [R.conditionID];

straight_idx = mod(conditionIDs, 3) == 1;
curved_idx   = mod(conditionIDs, 3) == 2;

straight_trials = find(straight_idx);
curved_trials   = find(curved_idx);
all_trials = [straight_trials(:); curved_trials(:)];
n_trials = numel(all_trials);
n_units = numel(R(1).unit);

windows.delay    = [-500, 0];
windows.movement = [0, 500];

window_labels = {'Delay Period', 'Movement Period'};
window_names  = {'delay', 'movement'};

n_windows = numel(window_labels);

bin_size = 50;
onset_align = [R.moveOnsetTime];

time_bins = windows.delay(1):bin_size:windows.movement(2);

bin_centers = time_bins(1:end-1) + bin_size/2;
n_bins = numel(bin_centers);

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

delay_bins = bin_centers >= -500 & bin_centers < 0;
movement_bins = bin_centers >= 0 & bin_centers < 500;

methods = {'svm', 'logreg'};
method_labels = {'SVM (Linear)', 'Logistic Regression'};

n_methods = numel(methods);
targetXY = reshape([R.targetXY], 2, [])';
angles = atan2d(targetXY(:,2), targetXY(:,1));
angleGroups = round(angles / 45) * 45;
angleGroups(angleGroups == 180) = -180;
uniqueGroups = unique(angleGroups);
nGroups = numel(uniqueGroups);

dir_labels_all = nan(numel(R),1);
for g = 1:nGroups
    dir_labels_all(angleGroups == uniqueGroups(g)) = g;
end

chance_level_dir = 1 / nGroups;
straight_pos = 1:numel(straight_trials);
curved_pos   = numel(straight_trials)+1 : n_trials;
barrier_conditions = { ...
    straight_pos, ...
    curved_pos};

for b = 1:n_barrier % sanity check to see how many trials per condition
    trial_ids = barrier_conditions{b};
    dir_labels = dir_labels_all(all_trials(trial_ids));
    fprintf('\n%s:\n', barrier_labels_str{b});
    tabulate(dir_labels)
end

barrier_labels_str = { ...
    'Straight (No Barriers)', ...
    'Curved (Barriers)'};

n_barrier = numel(barrier_conditions);
n_perm = 1000;
n_cv_folds = 5;
obs_acc_dir  = nan(n_methods, n_windows, n_barrier);
perm_acc_dir = nan(n_methods, n_windows, n_barrier, n_perm);
p_values_dir = nan(n_methods, n_windows, n_barrier);
pred_store = cell(n_methods, n_windows, n_barrier);
true_store = cell(n_methods, n_windows, n_barrier);

% ---- SANITY CHECK: all-trials decoding ----
X_all_delay    = squeeze(mean(firing_rates(:,:,delay_bins),    3));
X_all_movement = squeeze(mean(firing_rates(:,:,movement_bins), 3));
dir_labels_all_subset = dir_labels_all(all_trials);

cv_all = cvpartition(dir_labels_all_subset, 'KFold', 5, 'Stratify', true);

acc_all_delay    = classify_kfold(X_all_delay,    dir_labels_all_subset, 'svm', cv_all);
acc_all_movement = classify_kfold(X_all_movement, dir_labels_all_subset, 'svm', cv_all);

fprintf('\nAll-trials SVM decoding (delay):    %.3f\n', acc_all_delay);
fprintf('All-trials SVM decoding (movement): %.3f\n', acc_all_movement);
fprintf('Chance level: %.3f\n\n', 1/nGroups);
% -------------------------------------------

for b = 1:n_barrier
    trial_pos  = barrier_conditions{b};
    trial_ids  = all_trials(trial_pos); 
    dir_labels = dir_labels_all(trial_ids); 
    n_tr = numel(trial_pos);
    X_delay_b    = squeeze(mean(firing_rates(trial_pos,:,delay_bins),    3));
    X_movement_b = squeeze(mean(firing_rates(trial_pos,:,movement_bins), 3));
    window_features_b = {X_delay_b, X_movement_b};
    cv = cvpartition(dir_labels, 'KFold', n_cv_folds,'Stratify', true);

    for w = 1:n_windows
        X = window_features_b{w};
        for m = 1:n_methods
            fprintf('\n========================================\n');
            fprintf('%s | %s | %s\n', ...
                method_labels{m}, ...
                window_labels{w}, ...
                barrier_labels_str{b});
            fprintf('========================================\n');
            [pred, true_lab] = classify_kfold_predictions(X, dir_labels, methods{m}, cv);
            pred_store{m,w,b} = pred;
            true_store{m,w,b} = true_lab;
            obs_acc_dir(m,w,b) = mean(pred == true_lab);
            fprintf('Observed accuracy = %.3f\n', ...
                obs_acc_dir(m,w,b));
            method_m = methods{m};
            perm_acc_tmp = nan(n_perm,1);

            parfor p = 1:n_perm
                disp(p)
                shuf_labels = dir_labels(randperm(n_tr));
                perm_acc_tmp(p) = classify_kfold( ...
                    X, ...
                    shuf_labels, ...
                    method_m, ...
                    cv);
            end
            perm_acc_dir(m,w,b,:) = perm_acc_tmp;
            p_values_dir(m,w,b) = ...
                (sum(perm_acc_tmp >= obs_acc_dir(m,w,b)) + 1) ...
                / (n_perm + 1);

            fprintf('p-value = %.4f\n', ...
                p_values_dir(m,w,b));

        end
    end
end

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
                p_str = sprintf('%.4f', p_val);

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

%% Plotting
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
            pred = pred_store{m,w,b};
            true_lab = true_store{m,w,b};
            C = zeros(nGroups, nGroups);

            for i = 1:nGroups
                for j = 1:nGroups
                    C(i,j) = sum( ...
                        true_lab == i & pred == j);
                end
            end

            C_norm = C ./ sum(C,2);
            imagesc(ax, C_norm);
            colormap(ax, 'sky');

            clim(ax, [0 1]);
            hold(ax, 'on');
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
            for i = 1:nGroups
                rectangle(ax, ...
                    'Position', [i-0.5, i-0.5, 1, 1], ...
                    'EdgeColor', [0.9 0.5 0.1], ...
                    'LineWidth', 1.5);

            end
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

%% Helper functions
function acc = classify_kfold(X, y, method, cv)
    correct = 0;
    n_folds = cv.NumTestSets;
    for fold = 1:n_folds
        train_idx = training(cv, fold);
        test_idx = test(cv, fold);

        X_train = X(train_idx,:);
        y_train = y(train_idx);

        X_test = X(test_idx,:);
        y_test = y(test_idx);

        mu = mean(X_train,1);
        sigma = std(X_train,0,1);
        sigma(sigma == 0) = 1;
        X_train = (X_train - mu) ./ sigma;
        X_test = (X_test - mu) ./ sigma;

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
        pred = predict(mdl, X_test);
        correct = correct + sum(pred == y_test);
    end
    acc = correct / numel(y);
end

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
        mu = mean(X_train,1);
        sigma = std(X_train,0,1);
        sigma(sigma == 0) = 1;
        X_train = (X_train - mu) ./ sigma;
        X_test = (X_test - mu) ./ sigma;
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
        pred_all(test_idx) = predict(mdl, X_test);
    end
end

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