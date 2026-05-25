% Behavioral: Hand velocity over time split by the curved and straight 
% reaches. Compute peak hand speed during trial, split trials into groups 
% based on curved or straight reaches, compare the means. Monte Carlo: 
% shuffle the curved/straight labels, recompute the difference in the 
% means, do this 1,000 times to create a null distribution. P-value would 
% show us if curved and straight reaches have significantly different peak 
% hand speeds. 

%% Find peak speeds
function handVelocityMonteCarlo(R)
% handVelocityMonteCarlo  Compute peak hand speeds per trial, plot their
%   distributions for straight vs. curved conditions, and run a Monte Carlo
%   permutation test on the difference in means.
%
%   handVelocityMonteCarlo(R)
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
% Compute peak speed per trial
% -------------------------------------------------------------------------
peak_speeds = nan(1, num_trials);

for trial = 1:num_trials
    hand_times = R(trial).HAND.times(:);
    x_pos      = R(trial).HAND.X(:);
    y_pos      = R(trial).HAND.Y(:);

    [~, onset_idx] = min(abs(hand_times - moveOnsetTimes(trial)));
    [~, end_idx]   = min(abs(hand_times - moveEndsTimes(trial)));

    x_seg = x_pos(onset_idx:end_idx);
    y_seg = y_pos(onset_idx:end_idx);
    dt    = diff(hand_times(onset_idx:end_idx));
    dx    = diff(x_seg);
    dy    = diff(y_seg);

    valid = dt > 0;
    dt = dt(valid);
    dx = dx(valid);
    dy = dy(valid);

    speed = sqrt((dx ./ dt).^2 + (dy ./ dt).^2);
    peak_speeds(trial) = max(speed);
end

% -------------------------------------------------------------------------
% Split by condition, drop NaNs
% -------------------------------------------------------------------------
straight_speed = peak_speeds(straight_idx & ~isnan(peak_speeds));
curved_speed   = peak_speeds(curved_idx   & ~isnan(peak_speeds));

observed_diff = mean(curved_speed) - mean(straight_speed);

% -------------------------------------------------------------------------
% Figure 1: Speed distributions
% -------------------------------------------------------------------------
figure;

subplot(1, 2, 1);
histogram(straight_speed, 30, 'FaceColor', [0 0 0.8], ...
          'EdgeColor', 'none', 'Normalization', 'probability');
hold on;
xline(mean(straight_speed), 'b--', 'LineWidth', 2);
xlabel('Peak Hand Speed');
ylabel('Proportion of Trials');
title('No Barriers');

subplot(1, 2, 2);
histogram(curved_speed, 30, 'FaceColor', 'red', ...
          'EdgeColor', 'none', 'Normalization', 'probability');
hold on;
xline(mean(curved_speed), 'r--', 'LineWidth', 2);
xlabel('Peak Hand Speed');
ylabel('Proportion of Trials');
title('Barriers');

sgtitle('Peak Hand Speed Distribution by Condition');

% -------------------------------------------------------------------------
% Monte Carlo permutation test
% -------------------------------------------------------------------------
n_perm           = 10000;
permutation_diffs = nan(1, n_perm);
pooled_speeds    = [straight_speed, curved_speed];
num_straight     = numel(straight_speed);

for i = 1:n_perm
    shuffled     = pooled_speeds(randperm(numel(pooled_speeds)));
    perm_straight = shuffled(1:num_straight);
    perm_curved   = shuffled(num_straight+1:end);
    permutation_diffs(i) = mean(perm_curved) - mean(perm_straight);
end

p_value = mean(abs(permutation_diffs) >= abs(observed_diff));

% -------------------------------------------------------------------------
% Figure 2: Null distribution
% -------------------------------------------------------------------------
figure;
histogram(permutation_diffs, 50, 'FaceColor', [0.5 0.5 0.5], ...
          'EdgeColor', 'none', 'Normalization', 'probability');
hold on;
xline( observed_diff, 'r--', 'LineWidth', 2);
xline(-observed_diff, 'r--', 'LineWidth', 2);

xlabel('Difference in Mean Peak Speed (Curved - Straight)');
ylabel('Proportion of Permutations');
title(sprintf('Null Distribution (Monte Carlo Permutation Test)\nObserved difference = %.3f, p = %.3f', ...
              observed_diff, p_value));
legend('Null distribution', 'Observed difference');

end % handVelocityMonteCarlo