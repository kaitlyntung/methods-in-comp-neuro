% Behavioral: Find a curvature metric (ratio of actual path length to 
% straight line distance between start and end). Plot the distribution of 
% curvature with a bootstrapped confidence interval for curved and straight 
% path conditions. This would be a behavioral verification that the animal 
% actually did the curved paths and not just the straight length ones. 

function curvatureMetric(R)
% curvatureMetric  Compute path curvature ratios per trial, plot their
%   distributions for straight vs. curved conditions, and run a bootstrap
%   analysis on the difference in means with 95% CIs.
%
%   curvatureMetric(R)
%
%   Input:
%     R  - Trial data struct array. Each element R(tr) must contain:
%            .conditionID    - Scalar condition identifier
%            .moveOnsetTime  - Movement onset timestamp (ms)
%            .moveEndsTime   - Movement end timestamp (ms)
%            .HAND.times     - Hand position sample timestamps (ms)
%            .HAND.X         - Hand X positions
%            .HAND.Y         - Hand Y positions

% -------------------------------------------------------------------------
% Condition indexing
% -------------------------------------------------------------------------
num_trials   = numel(R);
conditionIDs = [R.conditionID];
straight_idx = mod(conditionIDs, 3) == 1;
curved_idx   = mod(conditionIDs, 3) == 2;

moveOnsetTimes = [R.moveOnsetTime];
moveEndsTimes  = [R.moveEndsTime];

% -------------------------------------------------------------------------
% Compute curvature ratio per trial
% -------------------------------------------------------------------------
curvature = nan(1, num_trials);

for trial = 1:num_trials
    hand_times = R(trial).HAND.times;
    x_pos      = R(trial).HAND.X;
    y_pos      = R(trial).HAND.Y;

    [~, onset_idx] = min(abs(hand_times - moveOnsetTimes(trial)));
    [~, end_idx]   = min(abs(hand_times - moveEndsTimes(trial)));

    x_seg = x_pos(onset_idx:end_idx);
    y_seg = y_pos(onset_idx:end_idx);

    coords            = [x_seg(:), y_seg(:)];
    segment_lengths   = sqrt(sum(diff(coords).^2, 2));
    total_path_length = sum(segment_lengths);
    straight_line_dist = sqrt((x_seg(end) - x_seg(1))^2 + ...
                               (y_seg(end) - y_seg(1))^2);

    curvature(trial) = total_path_length / straight_line_dist;
end

% -------------------------------------------------------------------------
% Split by condition, drop NaNs
% -------------------------------------------------------------------------
straight_curv = curvature(straight_idx & ~isnan(curvature));
curved_curv   = curvature(curved_idx   & ~isnan(curvature));

% -------------------------------------------------------------------------
% Figure 1: Curvature distributions
% -------------------------------------------------------------------------
figure;

subplot(1, 2, 1);
histogram(straight_curv, 30, 'FaceColor', [0 0 0.8], ...
          'EdgeColor', 'none', 'Normalization', 'probability');
xlabel('Curvature Ratio');
ylabel('Proportion of Trials');
title('No Barriers');
xlim([0.5 3]);

subplot(1, 2, 2);
histogram(curved_curv, 30, 'FaceColor', 'red', ...
          'EdgeColor', 'none', 'Normalization', 'probability');
xlabel('Curvature Ratio');
ylabel('Proportion of Trials');
title('Barriers');
xlim([0.5 3]);

sgtitle('Path Curvature Distribution by Condition');

% -------------------------------------------------------------------------
% Bootstrap
% -------------------------------------------------------------------------
n_bootstrap = 10000;

straight_random_inds = randi(numel(straight_curv), numel(straight_curv), n_bootstrap);
straight_bootstrap   = mean(straight_curv(straight_random_inds), 1);

curved_random_inds = randi(numel(curved_curv), numel(curved_curv), n_bootstrap);
curved_bootstrap   = mean(curved_curv(curved_random_inds), 1);

straight_mean = mean(straight_curv);
straight_ci   = prctile(straight_bootstrap, [2.5, 97.5]);

curved_mean = mean(curved_curv);
curved_ci   = prctile(curved_bootstrap, [2.5, 97.5]);

% -------------------------------------------------------------------------
% Figure 2: Bootstrap distributions with 95% CIs
% -------------------------------------------------------------------------
figure;
hold on;

histogram(straight_bootstrap, 50, 'FaceColor', [0 0 0.8], 'FaceAlpha', 0.5, ...
          'EdgeColor', 'none', 'Normalization', 'probability');
histogram(curved_bootstrap,   50, 'FaceColor', 'red',     'FaceAlpha', 0.5, ...
          'EdgeColor', 'none', 'Normalization', 'probability');

yl = ylim;

fill([straight_ci(1) straight_ci(2) straight_ci(2) straight_ci(1)], ...
     [0 0 yl(2) yl(2)], [0 0 0.8], 'FaceAlpha', 0.15, 'EdgeColor', 'none');
fill([curved_ci(1)   curved_ci(2)   curved_ci(2)   curved_ci(1)], ...
     [0 0 yl(2) yl(2)], 'red',     'FaceAlpha', 0.15, 'EdgeColor', 'none');

xline(straight_mean, 'b-', 'LineWidth', 2);
xline(curved_mean,   'r-', 'LineWidth', 2);

legend('No Barriers bootstrap', 'Barriers bootstrap', ...
       'No Barriers 95% CI',    'Barriers 95% CI', ...
       'No Barriers mean',      'Barriers mean', ...
       'Location', 'best');
xlabel('Mean Curvature Ratio');
ylabel('Proportion of Bootstrap Samples');
title('Bootstrap Mean Distributions with 95% CIs');

end % curvatureMetric