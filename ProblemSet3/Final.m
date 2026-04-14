%% Problem 1
colors = hsv(13);
figure; hold on;

for direction = 1:13
    for trial = 1:184
        spikes = find(mtNeuron.data(:,direction,trial));
        y = trial_offset + trial * ones(size(spikes));
        plot(spikes*2, y, '.', 'Color', colors(direction,:));
    end
    trial_offset = trial_offset + 184;
end

xlabel('Time (ms)');
ylabel('Trial');
title('Raster Plot');

%% Problem 2
cumulative_spikes = cumsum(mtNeuron.data, 1);
num_directions = size(mtNeuron.data, 2);
num_trials = size(mtNeuron.data, 3);
num_time = size(mtNeuron.data, 1);
% I(T) = I(T_n, theta) = sum_n P_T(n)log_s(P_T(n|theta)/P_T(n))

fractions = [1/4, 1/3, 1/2, 2/3, 1]; % eq7
num_fracs = length(fractions);
MI_raw = zeros(num_time, num_fracs);

for fi = 1:num_fracs
    n_samp = floor(fractions(fi) * num_trials); 
    for t = 1:num_time

        counts = zeros(num_directions, n_samp);
        for d = 1:num_directions
            idx = randperm(num_trials, n_samp);
            counts(d, :) = cumulative_spikes(t, d, idx);
        end

        all_counts = counts(:);
        max_n = max(all_counts);
        count_vals = 0:max_n;
        num_vals = length(count_vals);

        % P(n|theta) eq2 numerator
        Pn_given_theta = zeros(num_directions, num_vals);
        for d = 1:num_directions
            for ni = 1:num_vals
                Pn_given_theta(d, ni) = mean(counts(d, :) == count_vals(ni));
            end
        end
        % P(n) = (1/num_directions) * sum_theta P(n|theta) eq2
        Pn = mean(Pn_given_theta, 1);
        
        % ── Mutual information  I = sum_n P(n) log2[ P(n|theta) / P(n) ]
        % Equivalently: I = H(n) - <H(n|theta)>_theta  (Eqs. 1 & 5–6)
        % Computed as weighted average over directions of KL divergences

        MI = 0;
        P_theta = 1 / num_directions;   % uniform prior

        for d = 1:num_directions
            Pn_d = Pn_given_theta(d, :);    % P(n | theta_d)
            % Sum only where both Pn_d > 0 and Pn > 0 (avoid log(0))
            mask = Pn_d > 0 & Pn > 0;
            MI = MI + P_theta * sum( Pn_d(mask) .* log2( Pn_d(mask) ./ Pn(mask) ) );
        end

        MI_raw(t, fi) = MI;
    end
end

% ── Bias correction via linear extrapolation  (Eq. 7 in paper) ────────────
% I_est(N) = I_true + a/N  =>  plot I_raw vs 1/N and take y-intercept
MI_corrected = zeros(num_time, 1);
inv_fracs    = 1 ./ fractions;           % x-axis: 1/fraction (proxy for 1/N)

for t = 1:num_time
    p = polyfit(inv_fracs, MI_raw(t, :), 1);   % linear fit
    MI_corrected(t) = p(2);                     % y-intercept = I_true
end

MI_corrected = max(MI_corrected, 0);   % information can't be negative

% ── Time axis: 2 ms bins, referenced to stimulus onset ─────────────────────
time_ms = (0:num_time-1) * 2;   % 0 to 510 ms

% ── Plot ───────────────────────────────────────────────────────────────────
figure;
plot(time_ms, MI_corrected, 'k-', 'LineWidth', 1.5);
xlabel('Time from stimulus onset (ms)');
ylabel('Mutual information (bits)');
title('Mutual information between cumulative spike count and motion direction');
xlim([0 512]);
ylim([0 max(MI_corrected)*1.15]);
xline(256, '--', 'Stimulus offset', 'LabelVerticalAlignment','bottom');