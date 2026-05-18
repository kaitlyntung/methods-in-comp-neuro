windows.delay    = [-500, 0];
windows.movement = [0, 500];
window_labels    = {'Delay Period', 'Movement Period'};
window_names     = {'delay', 'movement'};
n_windows        = numel(window_labels);

all_trials = [straight_trials(:); curved_trials(:)];
labels     = [zeros(numel(straight_trials), 1); ones(numel(curved_trials), 1)];
n_trials   = numel(all_trials);

n_perm     = 1000;
n_cv_folds = 3;

delay_bins    = bin_centers >= -500 & bin_centers < 0;
movement_bins = bin_centers >= 0    & bin_centers < 500;

X_delay    = squeeze(mean(firing_rates(:, :, delay_bins),    3));
X_movement = squeeze(mean(firing_rates(:, :, movement_bins), 3));

window_features = {X_delay, X_movement};

methods       = {'svm', 'logreg'};
method_labels = {'SVM (Linear)', 'Logistic Regression'};
n_methods     = numel(methods);

obs_acc  = nan(n_methods, n_windows);
perm_acc = nan(n_methods, n_windows, n_perm);
p_values = nan(n_methods, n_windows);

%% --- Straight vs Curved classifier (the original loop that was missing) ---
for w = 1:n_windows
    X = window_features{w};
    for m = 1:n_methods
        obs_acc(m, w) = classify_kfold(X, labels, methods{m}, n_cv_folds);

        method_m     = methods{m};
        perm_acc_tmp = nan(n_perm, 1);
        parfor p = 1:n_perm
            shuf_labels     = labels(randperm(n_trials));
            perm_acc_tmp(p) = classify_kfold(X, shuf_labels, method_m, n_cv_folds);
        end

        perm_acc(m, w, :) = perm_acc_tmp;
        p_values(m, w)    = (sum(perm_acc_tmp >= obs_acc(m, w)) + 1) / (n_perm + 1);

        fprintf('Done: %s | %s\n', method_labels{m}, window_labels{w});
    end
end

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

%% --- Direction classifier within each barrier condition ---
targetXY    = reshape([R.targetXY], 2, [])';
angles      = atan2d(targetXY(:,2), targetXY(:,1));
angleGroups = round(angles / 45) * 45;
angleGroups(angleGroups == 180) = -180;
uniqueGroups = unique(angleGroups);
nGroups      = numel(uniqueGroups);

dir_labels_all = nan(numel(R), 1);
for g = 1:nGroups
    dir_labels_all(angleGroups == uniqueGroups(g)) = g;
end

chance_level_dir = 1 / nGroups;

barrier_conditions = {straight_trials, curved_trials};
barrier_labels_str = {'Straight (No Barriers)', 'Curved (Barriers)'};
n_barrier          = numel(barrier_conditions);

obs_acc_dir  = nan(n_methods, n_windows, n_barrier);
perm_acc_dir = nan(n_methods, n_windows, n_barrier, n_perm);
p_values_dir = nan(n_methods, n_windows, n_barrier);

% Build a lookup from trial ID -> row in firing_rates
% firing_rates rows correspond to all_trials order, so we verify this explicitly
fr_trial_index = all_trials;  % row i of firing_rates = trial fr_trial_index(i)

for b = 1:n_barrier
    trial_ids  = barrier_conditions{b};
    dir_labels = dir_labels_all(trial_ids);
    n_tr       = numel(trial_ids);

    % Map trial IDs to firing_rates rows safely
    [found, row_idx] = ismember(trial_ids, fr_trial_index);
    if any(~found)
        error('Some trials in barrier condition %d not found in firing_rates index.', b);
    end

    X_delay_b    = squeeze(mean(firing_rates(row_idx, :, delay_bins),    3));
    X_movement_b = squeeze(mean(firing_rates(row_idx, :, movement_bins), 3));
    window_features_b = {X_delay_b, X_movement_b};

    for w = 1:n_windows
        X = window_features_b{w};
        for m = 1:n_methods
            obs_acc_dir(m, w, b) = classify_kfold(X, dir_labels, methods{m}, n_cv_folds);

            method_m     = methods{m};
            perm_acc_tmp = nan(n_perm, 1);
            parfor p = 1:n_perm
                shuf_labels     = dir_labels(randperm(n_tr));
                perm_acc_tmp(p) = classify_kfold(X, shuf_labels, method_m, n_cv_folds);
            end

            perm_acc_dir(m, w, b, :) = perm_acc_tmp;
            p_values_dir(m, w, b)    = (sum(perm_acc_tmp >= obs_acc_dir(m,w,b)) + 1) / (n_perm + 1);

            fprintf('Done: %s | %s | %s\n', method_labels{m}, window_labels{w}, barrier_labels_str{b});
        end
    end
end

fprintf('\n%s\n', repmat('-', 1, 75));
fprintf('%-22s %-20s %-24s %-10s %-10s\n', 'Method', 'Window', 'Condition', 'Accuracy', 'p-value');
fprintf('%s\n', repmat('-', 1, 75));
for b = 1:n_barrier
    for w = 1:n_windows
        for m = 1:n_methods
            fprintf('%-22s %-20s %-24s %-10.3f %-10.4f %s\n', ...
                method_labels{m}, window_labels{w}, barrier_labels_str{b}, ...
                obs_acc_dir(m,w,b), p_values_dir(m,w,b), ...
                significance_label(p_values_dir(m,w,b)));
        end
    end
end
fprintf('%s\n', repmat('-', 1, 75));

%% Helper functions
function acc = classify_kfold(X, y, method, k)
    cv      = cvpartition(y, 'KFold', k, 'Stratify', true);
    correct = 0;
    for fold = 1:k
        X_train = X(cv.training(fold), :);
        y_train = y(cv.training(fold));
        X_test  = X(cv.test(fold),     :);
        y_test  = y(cv.test(fold));

        mu              = mean(X_train, 1);
        sigma           = std(X_train,  0, 1);
        sigma(sigma==0) = 1;
        X_train         = (X_train - mu) ./ sigma;
        X_test          = (X_test  - mu) ./ sigma;

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