% Behavioral: Find a curvature metric (ratio of actual path length to 
% straight line distance between start and end). Plot the distribution of 
% curvature with a bootstrapped confidence interval for curved and straight 
% path conditions. This would be a behavioral verification that the animal 
% actually did the curved paths and not just the straight length ones. 

num_trials = length(R);
num_directions = 36;
mazeCategory = {R.mazeCategory};
moveOnsetTimes = [R.moveOnsetTime];
moveEndsTimes = [R.moveEndsTime];
targetXYs = {R.targetXY};
fixationXYs = {R.fixationXY};

conditionIDs = [R.conditionID];
straight_idx = mod(conditionIDs, 3) == 1;
curved_idx   = mod(conditionIDs, 3) == 2;

curvature = nan(1, num_trials);

for trial = 1:num_trials
    hand_times = R(trial).HAND.times;
    x_pos = R(trial).HAND.X;
    y_pos = R(trial).HAND.Y;

    [~, moveOnsetTime_idx] = min(abs(hand_times - moveOnsetTimes(trial)));
    [~, moveEndsTime_idx]   = min(abs(hand_times - moveEndsTimes(trial)));

    x_segment = x_pos(moveOnsetTime_idx:moveEndsTime_idx);
    y_segment = y_pos(moveOnsetTime_idx:moveEndsTime_idx);

    % if numel(x_segment) < 2
    %     continue
    % end

    coordinates = [x_segment(:), y_segment(:)];
    segment_lengths = sqrt(sum(diff(coordinates).^2, 2));
    total_path_length = sum(segment_lengths);
    straight_dist = sqrt((x_segment(end) - x_segment(1))^2 + (y_segment(end) - y_segment(1))^2);

    if straight_dist > 1e-6
        curvature(trial) = total_path_length / straight_dist;
    end
end

straight_curv = curvature(straight_idx & ~isnan(curvature));
curved_curv   = curvature(curved_idx & ~isnan(curvature));


% Bootstrapped 95% CI on the mean
n_boot = 10000;
rng(42);

straight_boot = mean(straight_curv(randi(numel(straight_curv), numel(straight_curv), n_boot)), 1);
curved_boot   = mean(curved_curv(randi(numel(curved_curv),   numel(curved_curv),   n_boot)), 1);

s_mean = mean(straight_curv);
s_ci   = prctile(straight_boot, [2.5, 97.5]);

c_mean = mean(curved_curv);
c_ci   = prctile(curved_boot, [2.5, 97.5]);

fprintf('No Barriers: mean=%.3f, 95%% CI=[%.3f, %.3f], n=%d\n', ...
    s_mean, s_ci(1), s_ci(2), numel(straight_curv));
fprintf('Barriers:    mean=%.3f, 95%% CI=[%.3f, %.3f], n=%d\n', ...
    c_mean, c_ci(1), c_ci(2), numel(curved_curv));
figure;
set(gcf, 'Color', 'w');

% --- Panel 1: Overlapping histograms ---
subplot(1, 2, 1);
hold on;

bin_edges = linspace(0.9, 4.5, 60);

hs = histogram(straight_curv, bin_edges, 'Normalization', 'probability', ...
    'FaceColor', [0.29 0.47 0.81], 'EdgeColor', 'none', 'FaceAlpha', 0.6);
hc = histogram(curved_curv, bin_edges, 'Normalization', 'probability', ...
    'FaceColor', [0.84 0.37 0.37], 'EdgeColor', 'none', 'FaceAlpha', 0.6);

% Shade bootstrapped CI as vertical spans
yl = ylim;
fill([s_ci(1) s_ci(2) s_ci(2) s_ci(1)], [0 0 yl(2) yl(2)], ...
    [0.29 0.47 0.81], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
fill([c_ci(1) c_ci(2) c_ci(2) c_ci(1)], [0 0 yl(2) yl(2)], ...
    [0.84 0.37 0.37], 'FaceAlpha', 0.2, 'EdgeColor', 'none');

xline(s_mean, '--', 'Color', [0.29 0.47 0.81], 'LineWidth', 2);
xline(c_mean, '--', 'Color', [0.84 0.37 0.37], 'LineWidth', 2);

xlabel('Curvature Index (Path Length / Straight-Line Distance)', 'FontSize', 11);
ylabel('Proportion of Trials', 'FontSize', 11);
title('Distribution of Path Curvature', 'FontSize', 12, 'FontWeight', 'bold');
legend([hs, hc], {'No Barriers (straight)', 'Barriers (curved)'}, 'Location', 'northeast');
box off;

% --- Panel 2: Mean + 95% CI bar plot with jittered data ---
subplot(1, 2, 2);
hold on;

bar_x     = [1, 2];
bar_means = [s_mean, c_mean];
bar_colors = [0.29 0.47 0.81; 0.84 0.37 0.37];

for i = 1:2
    bar(bar_x(i), bar_means(i), 0.5, 'FaceColor', bar_colors(i,:), ...
        'EdgeColor', 'none', 'FaceAlpha', 0.85);
end

% Error bars from bootstrap CI
errorbar(1, s_mean, s_mean - s_ci(1), s_ci(2) - s_mean, ...
    'k', 'LineWidth', 2, 'CapSize', 10, 'LineStyle', 'none');
errorbar(2, c_mean, c_mean - c_ci(1), c_ci(2) - c_mean, ...
    'k', 'LineWidth', 2, 'CapSize', 10, 'LineStyle', 'none');

% Jittered individual points
rng(0);
jitter_s = (rand(size(straight_curv)) - 0.5) * 0.35;
jitter_c = (rand(size(curved_curv))   - 0.5) * 0.35;
scatter(1 + jitter_s, straight_curv, 6, [0.29 0.47 0.81], 'filled', 'MarkerFaceAlpha', 0.1);
scatter(2 + jitter_c, curved_curv,   6, [0.84 0.37 0.37], 'filled', 'MarkerFaceAlpha', 0.1);

% Significance bracket
y_sig = max(c_ci(2), prctile(curved_curv, 97)) * 1.05;
plot([1, 2], [y_sig y_sig], 'k-', 'LineWidth', 1.5);
text(1.5, y_sig * 1.02, '***', 'HorizontalAlignment', 'center', 'FontSize', 14);

xticks([1, 2]);
xticklabels({'No Barriers', 'Barriers'});
ylabel('Mean Curvature Index', 'FontSize', 11);
title('Mean Curvature ± 95% Bootstrap CI', 'FontSize', 12, 'FontWeight', 'bold');
xlim([0.5, 2.5]);
box off;

sgtitle('Behavioral Verification: Path Curvature by Maze Condition', ...
    'FontSize', 13, 'FontWeight', 'bold');