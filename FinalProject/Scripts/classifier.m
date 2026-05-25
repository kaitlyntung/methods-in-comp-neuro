% We will classify curved vs. straight reaches based on neural data from 
% the delay and movement period separately using SVM and logistic 
% regression to assess when kinematic information is encoded in M1.

% We will assess the significance of our classifier using a Monte Carlo 
% shuffle test, where we shuffle curved vs straight trial labels 1000 times
% to generate a null distribution and compute a p-value for each
% method/time-window combination.

function results = classifier(R, varargin)

p = inputParser;
addParameter(p, 'BinSize',  50,              @isnumeric);
addParameter(p, 'DelayWin', [-500  0],       @(x) isnumeric(x) && numel(x)==2);
addParameter(p, 'MoveWin',  [0    500],      @(x) isnumeric(x) && numel(x)==2);
addParameter(p, 'nPerm',    1000,            @isnumeric);
addParameter(p, 'nFolds',   5,               @isnumeric);
addParameter(p, 'Methods',  {'svm','logreg'},@iscell);
parse(p, varargin{:});
opt = p.Results;

bin_size     = opt.BinSize;
delay_win    = opt.DelayWin;
move_win     = opt.MoveWin;
n_perm       = opt.nPerm;
n_cv_folds   = opt.nFolds;
methods      = opt.Methods;
method_labels = cellfun(@(m) method_display_name(m), methods, ...
                        'UniformOutput', false);
n_methods    = numel(methods);

window_labels = {'Delay Period', 'Movement Period'};
n_windows     = numel(window_labels);

barrier_labels_str = {'Straight (No Barriers)', 'Curved (Barriers)'};
n_barrier          = numel(barrier_labels_str);

angle_strs = {'Left','Lower-Left','Lower-Right','Upper-Right','Up','Upper-Left'};

% ── Index straight vs curved trials ──────────────────────────────────────
conditionIDs  = [R.conditionID];
straight_trials = find(mod(conditionIDs, 3) == 1);
curved_trials   = find(mod(conditionIDs, 3) == 2);
all_trials      = [straight_trials(:); curved_trials(:)];
n_trials        = numel(all_trials);
n_units         = numel(R(1).unit);

time_bins   = delay_win(1) : bin_size : move_win(2);
bin_centers = time_bins(1:end-1) + bin_size/2;
n_bins      = numel(bin_centers);

onset_align  = [R.moveOnsetTime];
firing_rates = nan(n_trials, n_units, n_bins);

for t = 1:n_trials
    tr = all_trials(t);
    t0 = onset_align(tr);
    for u = 1:n_units
        spike_times = R(tr).unit(u).spikeTimes - t0;
        counts      = histcounts(spike_times, time_bins);
        firing_rates(t, u, :) = counts / (bin_size / 1000);
    end
end
disp('Firing rates computed.');

delay_bins    = bin_centers >= delay_win(1) & bin_centers < delay_win(2);
movement_bins = bin_centers >= move_win(1)  & bin_centers < move_win(2);

targetXY    = reshape([R.targetXY], 2, [])';
angles      = atan2d(targetXY(:,2), targetXY(:,1));
angleGroups = round(angles / 45) * 45;
angleGroups(angleGroups ==  180)             = -180;
angleGroups(angleGroups == 0 & angles <  5)  = -45;
angleGroups(angleGroups == 0 & angles >= 5)  =  45;
uniqueGroups = unique(angleGroups);
nGroups      = length(uniqueGroups);

dir_labels_all = nan(numel(R), 1);
for g = 1:nGroups
    dir_labels_all(angleGroups == uniqueGroups(g)) = g;
end
chance_level_dir = 1 / nGroups;

straight_pos = 1 : numel(straight_trials);
curved_pos   = numel(straight_trials)+1 : n_trials;
barrier_conditions = {straight_pos, curved_pos};

for b = 1:n_barrier
    dir_labels = dir_labels_all(all_trials(barrier_conditions{b}));
    fprintf('\n%s:\n', barrier_labels_str{b});
    tabulate(dir_labels);
end

obs_acc  = nan(n_methods, n_windows, n_barrier);
perm_acc = nan(n_methods, n_windows, n_barrier, n_perm);
p_values = nan(n_methods, n_windows, n_barrier);
pred_store = cell(n_methods, n_windows, n_barrier);
true_store = cell(n_methods, n_windows, n_barrier);

for b = 1:n_barrier
    trial_pos  = barrier_conditions{b};
    trial_ids  = all_trials(trial_pos); 
    dir_labels = dir_labels_all(trial_ids);
    n_tr       = numel(trial_pos);

    X_delay    = squeeze(mean(firing_rates(trial_pos, :, delay_bins),    3));
    X_movement = squeeze(mean(firing_rates(trial_pos, :, movement_bins), 3));
    window_features = {X_delay, X_movement};

    cv = cvpartition(dir_labels, 'KFold', n_cv_folds, 'Stratify', true);

    for w = 1:n_windows
        X = window_features{w};

        for m = 1:n_methods
            fprintf('\n========================================\n');
            fprintf('%s | %s | %s\n', method_labels{m}, ...
                window_labels{w}, barrier_labels_str{b});
            fprintf('========================================\n');

            [pred, true_lab] = dd_classify_kfold_predictions( ...
                X, dir_labels, methods{m}, cv);

            pred_store{m,w,b} = pred;
            true_store{m,w,b} = true_lab;
            obs_acc(m,w,b)    = mean(pred == true_lab);
            fprintf('Observed accuracy = %.3f\n', obs_acc(m,w,b));

            % Permutation test
            method_m      = methods{m};
            perm_acc_tmp  = nan(n_perm, 1);
            parfor pp = 1:n_perm
                shuf = dir_labels(randperm(n_tr));
                perm_acc_tmp(pp) = dd_classify_kfold( ...
                    X, shuf, method_m, cv);
            end
            perm_acc(m,w,b,:) = perm_acc_tmp;
            p_values(m,w,b) = ...
                (sum(perm_acc_tmp >= obs_acc(m,w,b)) + 1) / (n_perm + 1);
            fprintf('p-value = %.4f\n', p_values(m,w,b));
        end
    end
end

fprintf('\n%s\n', repmat('-',1,75));
fprintf('%-22s %-20s %-24s %-10s %-10s\n', ...
    'Method','Window','Condition','Accuracy','p-value');
fprintf('%s\n', repmat('-',1,75));
for b = 1:n_barrier
    for w = 1:n_windows
        for m = 1:n_methods
            pv = p_values(m,w,b);
            fprintf('%-22s %-20s %-24s %-10.3f %-10s %s\n', ...
                method_labels{m}, window_labels{w}, barrier_labels_str{b}, ...
                obs_acc(m,w,b), sprintf('%.4f',pv), dd_significance_label(pv));
        end
    end
end
fprintf('%s\n', repmat('-',1,75));

plot_confusion_matrices(pred_store, true_store, p_values, obs_acc, ...
    n_barrier, n_methods, n_windows, nGroups, ...
    barrier_labels_str, method_labels, window_labels, angle_strs);

results.obs_acc   = obs_acc;
results.perm_acc  = perm_acc;
results.p_values  = p_values;
results.pred      = pred_store;
results.true_lab  = true_store;
results.params    = opt;

end 

function plot_confusion_matrices(pred_store, true_store, p_values, obs_acc, ...
        n_barrier, n_methods, n_windows, nGroups, ...
        barrier_labels_str, method_labels, window_labels, angle_strs)

    for b = 1:n_barrier
        figure('Color','w', ...
               'Position',[100 100 1200 1000], ...
               'Name', barrier_labels_str{b});
        panel = 0;

        for m = 1:n_methods
            for w = 1:n_windows
                panel = panel + 1;
                ax    = subplot(n_methods, n_windows, panel);

                pred     = pred_store{m,w,b};
                true_lab = true_store{m,w,b};

                % Build raw confusion matrix
                C = zeros(nGroups, nGroups);
                for i = 1:nGroups
                    for j = 1:nGroups
                        C(i,j) = sum(true_lab == i & pred == j);
                    end
                end
                C_norm = C ./ sum(C, 2);

                imagesc(ax, C_norm);
                colormap(ax, 'sky');
                clim(ax, [0 1]);
                hold(ax, 'on');

                % Cell text annotations
                for i = 1:nGroups
                    for j = 1:nGroups
                        val     = C_norm(i,j);
                        txt_col = 'k';
                        if val > 0.6; txt_col = 'w'; end
                        text(ax, j, i, sprintf('%.0f%%', val*100), ...
                            'HorizontalAlignment','center', ...
                            'VerticalAlignment','middle', ...
                            'FontSize', 8, 'Color', txt_col);
                    end
                end

                % Diagonal highlight boxes
                for i = 1:nGroups
                    rectangle(ax, ...
                        'Position',  [i-0.5, i-0.5, 1, 1], ...
                        'EdgeColor', [0.9 0.5 0.1], ...
                        'LineWidth', 1.5);
                end

                % Axis labels / title
                pv  = p_values(m,w,b);
                if pv < 0.0001
                    p_str = '< 0.0001';
                else
                    p_str = sprintf('%.4f', pv);
                end
                overall_acc = obs_acc(m,w,b) * 100;

                ax.XTick           = 1:nGroups;
                ax.YTick           = 1:nGroups;
                ax.XTickLabel      = angle_strs;
                ax.YTickLabel      = angle_strs;
                ax.TickLength      = [0 0];
                ax.FontSize        = 9;
                ax.XTickLabelRotation = 45;
                xlabel(ax, 'Predicted direction');
                ylabel(ax, 'True direction');
                title(ax, sprintf('%s | %s\nacc = %.1f%%, p = %s %s', ...
                    method_labels{m}, window_labels{w}, ...
                    overall_acc, p_str, dd_significance_label(pv)));

                cb               = colorbar(ax);
                cb.Label.String  = 'Proportion';
                cb.FontSize      = 8;
            end
        end

        sgtitle(sprintf('Direction Confusion Matrix — %s', ...
            barrier_labels_str{b}), ...
            'FontSize', 13, 'FontWeight', 'bold');
    end
end


function acc = dd_classify_kfold(X, y, method, cv)
% DD_CLASSIFY_KFOLD  K-fold accuracy for permutation testing.
    correct = 0;
    for fold = 1:cv.NumTestSets
        [X_tr, X_te, y_tr, y_te] = split_fold(X, y, cv, fold);
        mdl     = fit_model(X_tr, y_tr, method);
        pred    = predict(mdl, X_te);
        correct = correct + sum(pred == y_te);
    end
    acc = correct / numel(y);
end


function [pred_all, true_all] = dd_classify_kfold_predictions(X, y, method, cv)
% DD_CLASSIFY_KFOLD_PREDICTIONS  K-fold predictions (for confusion matrix).
    pred_all = nan(size(y));
    true_all = y;
    for fold = 1:cv.NumTestSets
        test_idx            = test(cv, fold);
        [X_tr, X_te, y_tr, ~] = split_fold(X, y, cv, fold);
        mdl                 = fit_model(X_tr, y_tr, method);
        pred_all(test_idx)  = predict(mdl, X_te);
    end
end


function [X_train, X_test, y_train, y_test] = split_fold(X, y, cv, fold)
% SPLIT_FOLD  Extract and z-score one CV fold.
    train_idx = training(cv, fold);
    test_idx  = test(cv,     fold);

    X_train = X(train_idx, :);
    X_test  = X(test_idx,  :);
    y_train = y(train_idx);
    y_test  = y(test_idx);

    mu              = mean(X_train, 1);
    sigma           = std(X_train,  0, 1);
    sigma(sigma==0) = 1;
    X_train = (X_train - mu) ./ sigma;
    X_test  = (X_test  - mu) ./ sigma;
end


function mdl = fit_model(X_train, y_train, method)
% FIT_MODEL  Fit a multiclass ECOC model.
    switch method
        case 'svm'
            t = templateSVM('KernelFunction', 'rbf', ...
                'Standardize',    false, ...
                'BoxConstraint',  1, ...
                'KernelScale',    'auto');
        case 'logreg'
            t = templateLinear('Learner',       'logistic', ...
                               'Regularization','ridge', ...
                               'Lambda',        1/size(X_train,1));
        otherwise
            error('Unknown method: %s', method);
    end
    mdl = fitcecoc(X_train, y_train, 'Learners', t);
end


function label = dd_significance_label(p)
% DD_SIGNIFICANCE_LABEL  Return significance asterisks.
    if     p < 0.001; label = '***';
    elseif p < 0.01;  label = '**';
    elseif p < 0.05;  label = '*';
    else;             label = '(n.s.)';
    end
end


function name = method_display_name(method)
% METHOD_DISPLAY_NAME  Human-readable label for each method string.
    switch method
        case 'svm';    name = 'SVM (RBF)';
        case 'logreg'; name = 'Logistic Regression';
        otherwise;     name = method;
    end
end