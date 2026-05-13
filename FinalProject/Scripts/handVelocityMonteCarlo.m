% Behavioral: Hand velocity over time split by the curved and straight 
% reaches. Compute peak hand speed during trial, split trials into groups 
% based on curved or straight reaches, compare the means. Monte Carlo: 
% shuffle the curved/straight labels, recompute the difference in the 
% means, do this 1,000 times to create a null distribution. P-value would 
% show us if curved and straight reaches have significantly different peak 
% hand speeds. 

peak_speed = nan(1, num_trials);

for trial = 1:num_trials
    hand_times = R(trial).HAND.times(:);
    x_pos      = R(trial).HAND.X(:);
    y_pos      = R(trial).HAND.Y(:);

    [~, onset_idx] = min(abs(hand_times - moveOnsetTimes(trial)));
    [~, end_idx]   = min(abs(hand_times - moveEndsTimes(trial)));

    x_seg = x_pos(onset_idx:end_idx);
    y_seg = y_pos(onset_idx:end_idx);

    if numel(x_seg) < 2
        continue
    end

    dt = diff(hand_times(onset_idx:end_idx));
    dx = diff(x_seg);
    dy = diff(y_seg);

    valid = dt > 0;
    dt = dt(valid);
    dx = dx(valid);
    dy = dy(valid);

    if isempty(dt)
        continue
    end

    speed = sqrt((dx ./ dt).^2 + (dy ./ dt).^2);
    peak_speed(trial) = max(speed);
end

% Split by condition
straight_speed = peak_speed(straight_idx & ~isnan(peak_speed));
curved_speed   = peak_speed(curved_idx   & ~isnan(peak_speed));

observed_diff = mean(curved_speed) - mean(straight_speed);

% Monte Carlo permutation test
n_perm     = 1000;
perm_diffs = nan(1, n_perm);
pooled     = [straight_speed, curved_speed];
n_straight = numel(straight_speed);

rng(42);
for i = 1:n_perm
    shuffled       = pooled(randperm(numel(pooled)));
    perm_straight  = shuffled(1:n_straight);
    perm_curved    = shuffled(n_straight+1:end);
    perm_diffs(i)  = mean(perm_curved) - mean(perm_straight);
end

p_value = mean(abs(perm_diffs) >= abs(observed_diff));

fprintf('Mean peak speed - No Barriers: %.2f  Barriers: %.2f\n', ...
    mean(straight_speed), mean(curved_speed));
fprintf('Observed difference (curved - straight): %.2f\n', observed_diff);
fprintf('Monte Carlo p-value (two-tailed, n=%d permutations): %.4f\n', n_perm, p_value);

% --- Plot ---
figure;
set(gcf, 'Color', 'w');

% Panel 1: Null distribution with observed difference
subplot(1, 2, 1);
hold on;

histogram(perm_diffs, 40, 'Normalization', 'probability', ...
    'FaceColor', [0.6 0.6 0.6], 'EdgeColor', 'none', 'FaceAlpha', 0.8);
xline(observed_diff,  'r-',  'LineWidth', 2.5, 'DisplayName', 'Observed diff');
xline(-observed_diff, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Mirrored (two-tailed)');

xlabel('Difference in Mean Peak Speed (curved - straight)', 'FontSize', 11);
ylabel('Proportion', 'FontSize', 11);
title('Monte Carlo Null Distribution', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'northeast');
text(0.05, 0.92, sprintf('p = %.4f', p_value), 'Units', 'normalized', ...
    'FontSize', 11, 'Color', 'r');
box off;

% Panel 2: Mean peak speed per condition with jittered points
subplot(1, 2, 2);
hold on;

bar_colors = [0.29 0.47 0.81; 0.84 0.37 0.37];
bar_means  = [mean(straight_speed), mean(curved_speed)];
bar_sems   = [std(straight_speed)/sqrt(numel(straight_speed)), ...
              std(curved_speed)/sqrt(numel(curved_speed))];

for i = 1:2
    bar(i, bar_means(i), 0.5, 'FaceColor', bar_colors(i,:), ...
        'EdgeColor', 'none', 'FaceAlpha', 0.85);
end

errorbar([1, 2], bar_means, bar_sems, 'k', 'LineWidth', 2, ...
    'CapSize', 10, 'LineStyle', 'none');

rng(0);
scatter(1 + (rand(size(straight_speed)) - 0.5) * 0.35, straight_speed, ...
    6, bar_colors(1,:), 'filled', 'MarkerFaceAlpha', 0.15);
scatter(2 + (rand(size(curved_speed))   - 0.5) * 0.35, curved_speed, ...
    6, bar_colors(2,:), 'filled', 'MarkerFaceAlpha', 0.15);

% Significance bracket
y_sig = max([bar_means + bar_sems]) * 1.1;
plot([1, 2], [y_sig y_sig], 'k-', 'LineWidth', 1.5);
if p_value < 0.001
    sig_str = '***';
elseif p_value < 0.01
    sig_str = '**';
elseif p_value < 0.05
    sig_str = '*';
else
    sig_str = 'n.s.';
end
text(1.5, y_sig * 1.02, sig_str, 'HorizontalAlignment', 'center', 'FontSize', 14);

xticks([1, 2]);
xticklabels({'No Barriers', 'Barriers'});
ylabel('Mean Peak Hand Speed', 'FontSize', 11);
title('Peak Hand Speed by Condition', 'FontSize', 12, 'FontWeight', 'bold');
xlim([0.5, 2.5]);
box off;

sgtitle('Hand Speed Analysis: Curved vs Straight Reaches', ...
    'FontSize', 13, 'FontWeight', 'bold');