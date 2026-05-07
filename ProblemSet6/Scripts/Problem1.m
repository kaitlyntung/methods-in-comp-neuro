%% Part A
start_time = 100;
end_time = 400;


num_trials = length(R);
num_units = length(R(1).unit);
spikeCounts = zeros(num_units, num_trials);

for trial = 1:num_trials
    targetAppearsTime = R(trial).targetAppearsTime;
    start_window = targetAppearsTime + start_time;
    end_window = targetAppearsTime + end_time;
    for unit = 1:num_units
        spike = R(trial).unit(unit).spikeTimes;
        spikeCounts(unit, trial) = sum(spike > start_window & spike <= end_window);
    end
end
%% Part B
window_size = (abs(end_time) + abs(start_time)) * 0.001;

conditions = [R.conditionID];

cond_one = conditions(1);
cond_two = conditions(2);

cond_one_idx = conditions == cond_one;
cond_two_idx = conditions == cond_two;

avgRateOne = mean(spikeCounts(:, cond_one_idx), 2) / window_size;
avgRateTwo = mean(spikeCounts(:, cond_two_idx), 2) / window_size;

%% Part C
diffRate = abs(avgRateOne - avgRateTwo);

tunedNeuronMask = diffRate >= 2;
tunedNeuronIdx = find(tunedNeuronMask);

fprintf('\nPart 1C — Tuned units (|ΔFR| ≥ 2 sp/s): %d / %d\n', ...
        sum(tunedNeuronMask), num_units);

fprintf('\n%5s  %14s  %14s  %8s\n', 'Unit', ...
        sprintf('Cond%d (sp/s)', cond_one), ...
        sprintf('Cond%d (sp/s)', cond_two), '|Diff|');

fprintf('%s\n', repmat('-', 1, 48));

for k = 1:length(tunedNeuronIdx)
    u = tunedNeuronIdx(k);

    fprintf('%5d  %14.2f  %14.2f  %8.2f\n', ...
            u, avgRateOne(u), avgRateTwo(u), diffRate(u));
end

%% Part 2
X = spikeCounts(tunedNeuronIdx, :)';
y = conditions';
 
mdl = fitcnb(X, y, 'DistributionNames', 'kernel');
 
cv = crossval(mdl, 'KFold', 10);
lossCV = kfoldLoss(cv); 
accCV = (1 - lossCV) * 100;
 
fprintf('\n── Part 2 Results ───────────────────────────────────\n');
fprintf('10-fold CV misclassification rate : %.4f\n', lossCV);
fprintf('10-fold CV accuracy               : %.2f%%\n', accCV);

%% Part 3A: Logistic Regression with glmfit + 10-fold CV
% Recode labels as 0 and 1 (required for binomial glmfit)
y_bin = double(conditions' == cond_two);   % cond_two → 1, cond_one → 0
 
X_tuned = spikeCounts(tunedNeuronIdx, :)'; % nTrials x nTunedUnits
 
% 10-fold cross-validation partition (fixed seed for reproducibility)
rng(42);
cv3 = cvpartition(y_bin, 'KFold', 10);
 
correctLR = 0;
for fold = 1:cv3.NumTestSets
    % Training and test sets for this fold
    X_train = X_tuned(cv3.training(fold), :);
    y_train = y_bin(cv3.training(fold));
    X_test  = X_tuned(cv3.test(fold), :);
    y_test  = y_bin(cv3.test(fold));
 
    % Fit logistic regression on training data
    % glmfit prepends an intercept column automatically
    b = glmfit(X_train, y_train, 'binomial', 'link', 'logit');
 
    % Predict probabilities on test data, then threshold at 0.5
    pHat = glmval(b, X_test, 'logit');
    yPred = double(pHat >= 0.5);
 
    correctLR = correctLR + sum(yPred == y_test);
end
 
accLR = correctLR / num_trials * 100;
 
fprintf('\n── Part 3A Results ──────────────────────────────────\n');
fprintf('10-fold CV accuracy (Logistic Regression): %.2f%%\n', accLR);
 
% ── Answer: Why does glmfit produce warnings? ─────────────────────────────────
% With 28 features and only ~52 training trials per fold, the data are likely
% perfectly (or near-perfectly) separable in the high-dimensional feature space.
% When this happens, the maximum-likelihood estimates of the logistic regression
% coefficients diverge to ±infinity — the algorithm can push the decision
% boundary to perfectly separate the training classes, so no finite set of
% coefficients truly maximises the likelihood. glmfit detects this as a failure
% of the IRLS algorithm to converge and warns accordingly. The model is
% overparameterised: too many free weights relative to the number of training
% examples.
 
%% Part 3B: ElasticNet-penalized Logistic Regression with lassoglm
% Alpha = 0.1  →  90% L2, 10% L1 (nearly ridge regression)
% Use built-in CV inside lassoglm to choose Lambda, then evaluate on a held-out
% test set that was never seen during fitting or lambda selection.
 
% Hold out 20% of trials as a final test set
rng(42);
cvOuter = cvpartition(y_bin, 'HoldOut', 0.2);
 
X_fit  = X_tuned(cvOuter.training, :);   % 80% for fitting + lambda selection
y_fit  = y_bin(cvOuter.training);
X_test = X_tuned(cvOuter.test, :);       % 20% held-out test set
y_test = y_bin(cvOuter.test);
 
% lassoglm with built-in 10-fold CV to select Lambda
% 'CV', 10  →  choose Lambda that minimises deviance across 10 inner folds
[B, FitInfo] = lassoglm(X_fit, y_fit, 'binomial', ...
                         'Alpha', 0.1, ...
                         'CV', 10);
 
% Retrieve the coefficients at the Lambda with minimum CV deviance
lambdaIdx   = FitInfo.IndexMinDeviance;
B_best      = B(:, lambdaIdx);
intercept   = FitInfo.Intercept(lambdaIdx);
 
% Predict on held-out test set
logits = X_test * B_best + intercept;
pHat   = 1 ./ (1 + exp(-logits));
yPred  = double(pHat >= 0.5);
 
accEN = sum(yPred == y_test) / length(y_test) * 100;
 
fprintf('\n── Part 3B Results ──────────────────────────────────\n');
fprintf('Best Lambda (min CV deviance) : %.4f\n', FitInfo.LambdaMinDeviance);
fprintf('Non-zero coefficients         : %d / %d\n', sum(B_best ~= 0), length(B_best));
fprintf('Test set accuracy (ElasticNet): %.2f%%\n', accEN);
 
% ── Answer: Why does regularization make the warnings go away? ────────────────
% The ElasticNet penalty adds a term to the cost function that shrinks
% coefficients toward zero (L2 component) and can zero some out entirely
% (L1 component). This prevents coefficients from diverging to ±infinity even
% when the data are perfectly separable, because any increase in coefficient
% magnitude is now directly penalised. The optimisation therefore always has a
% finite, well-defined solution and IRLS converges without warnings. The
% regularisation effectively constrains the model complexity so it cannot
% perfectly memorise the training data, trading a tiny amount of training-set
% fit for a stable, finite solution.
 
%% Part 4A: Linear SVM
% SVMs use +1 / -1 labels
y_svm = ones(num_trials, 1);
y_svm(conditions' == cond_one) = -1;   % cond_one → -1, cond_two → +1
 
mdlSVMlin = fitcsvm(X_tuned, y_svm, 'KernelFunction', 'linear');
 
cvSVMlin  = crossval(mdlSVMlin, 'KFold', 10);
accSVMlin = (1 - kfoldLoss(cvSVMlin)) * 100;
 
fprintf('\n── Part 4A Results ──────────────────────────────────\n');
fprintf('10-fold CV accuracy (Linear SVM): %.2f%%\n', accSVMlin);
 
% ── Answer: Why no convergence trouble? ───────────────────────────────────────
% Logistic regression finds the maximum-likelihood decision boundary, and when
% classes are perfectly separable the likelihood keeps increasing as coefficients
% grow toward ±infinity — there is no finite optimum, so IRLS never converges.
% A linear SVM is explicitly designed for separable data: instead of maximising
% the likelihood it maximises the margin between the two nearest points
% (the support vectors). Even with perfect separability this margin has a unique,
% finite maximum, so the quadratic program always has a well-defined solution and
% no warnings are produced.
 
%% Part 4B: RBF SVM with automatic hyperparameter optimisation
% OptimizeHyperparameters='auto' searches over BoxConstraint and KernelScale.
% HyperparameterOptimizationOptions sets the inner CV used during search.
rng(42);
mdlSVMrbf = fitcsvm(X_tuned, y_svm, ...
    'KernelFunction',              'rbf', ...
    'OptimizeHyperparameters',     'auto', ...
    'HyperparameterOptimizationOptions', ...
        struct('AcquisitionFunctionName', 'expected-improvement-plus', ...
               'ShowPlots', false, ...
               'Verbose',   0));
 
cvSVMrbf  = crossval(mdlSVMrbf, 'KFold', 10);
accSVMrbf = (1 - kfoldLoss(cvSVMrbf)) * 100;
 
fprintf('\n── Part 4B Results ──────────────────────────────────\n');
fprintf('10-fold CV accuracy (RBF SVM): %.2f%%\n', accSVMrbf);
 
% ── Answer: How does RBF SVM compare, and why? ────────────────────────────────
% The RBF kernel maps the data into an infinite-dimensional feature space,
% giving the SVM much more expressive power to fit non-linear boundaries.
% However, with only 58 trials and 28 features, that extra flexibility is more
% likely to overfit than to generalise — the model can memorise the training
% data but not necessarily the underlying structure. As a result the RBF SVM
% often performs similarly to or slightly worse than the linear SVM under CV,
% even after hyperparameter tuning. This illustrates the bias-variance tradeoff:
% more complex models do not always win when data are scarce.
 
%% Part 5: Fitting a Dynamical System
% Uses mazeJ0918Simple.mat (same R struct format as PS5).
% Parts A-D: epoch -150 to +200 ms around movement onset.
% Part E:    epoch -150 to  -50 ms around movement onset (preparatory).

%% ── Shared parameters ────────────────────────────────────────────────────────
SD_MS      = 30;       % Gaussian kernel SD for smoothing (ms)
BIN_MS     = 1;        % 1 ms bins
N_PCS      = 6;        % number of PCs to keep

% We'll run the full pipeline twice (Parts A-D, then Part E),
% so wrap it in a helper at the bottom of the file.

%% ── Parts A-D: Full movement epoch (-150 to +200 ms) ────────────────────────
epoch_AD = [-150, 200];
[M_AD, eigvals_AD, ~] = fitDynamics(R, epoch_AD, SD_MS, BIN_MS, N_PCS);

fprintf('── Parts A-D  (epoch %d to %+d ms) ──────────────────\n', epoch_AD);
fprintf('M matrix estimated (%dx%d)\n', size(M_AD,1), size(M_AD,2));
fprintf('Eigenvalues of M:\n');
disp(eigvals_AD);

%% ── Part E: Preparatory sub-epoch (-150 to -50 ms) ──────────────────────────
epoch_E = [-150, -50];
[M_E, eigvals_E, ~] = fitDynamics(R, epoch_E, SD_MS, BIN_MS, N_PCS);

fprintf('── Part E  (epoch %d to %+d ms) ──────────────────────\n', epoch_E);
fprintf('Eigenvalues of M:\n');
disp(eigvals_E);

%% ── Part D + E: Plot eigenvalues on the complex plane ────────────────────────
figure; hold on; grid on; axis equal;

% Unit circle for reference
theta = linspace(0, 2*pi, 300);
plot(cos(theta), sin(theta), 'k--', 'LineWidth', 0.8, 'DisplayName', 'Unit circle');

% Full epoch eigenvalues (blue)
scatter(real(eigvals_AD), imag(eigvals_AD), 80, 'b', 'filled', ...
        'DisplayName', sprintf('Full epoch (%d to %+d ms)', epoch_AD));

% Preparatory sub-epoch eigenvalues (red)
scatter(real(eigvals_E), imag(eigvals_E), 80, 'r', 'filled', ...
        'DisplayName', sprintf('Pre-movement (%d to %+d ms)', epoch_E));

xline(0, 'k:', 'LineWidth', 0.8); yline(0, 'k:', 'LineWidth', 0.8);
xlabel('Re(\lambda)');
ylabel('Im(\lambda)');
title('Eigenvalues of fitted linear dynamical system M');
legend('Location', 'best');

% ── Answer: How are the eigenvalues different? ────────────────────────────────
%
% The full-movement-epoch eigenvalues (blue) tend to have imaginary parts —
% they appear as conjugate pairs away from the real axis, indicating rotational
% (oscillatory) dynamics. This matches the rotational structure seen in motor
% cortex during movement: neural trajectories trace out loops in PC space,
% which a linear system reproduces through complex-eigenvalue rotation.
%
% The pre-movement (preparatory) eigenvalues (red) tend to lie closer to the
% real axis, with smaller imaginary parts, indicating dynamics that are more
% purely expanding or contracting rather than rotating. Preparatory activity
% is thought to set up an initial condition — the neural state is being
% "loaded" toward a particular starting point for the upcoming movement, a
% process that does not require oscillation. Once the Go cue arrives the
% system transitions into the rotational regime that drives the movement
% itself. This is consistent with the class figure showing that eigenvalues
% with large imaginary parts (rotation) appear specifically during the
% movement epoch, while the preparatory epoch is dominated by real eigenvalues
% (growth/decay along fixed axes).

%% ══════════════════════════════════════════════════════════════════════════════
%% Helper function: smooth → trial-average → PCA → estimate M via regression
%% ══════════════════════════════════════════════════════════════════════════════
function [M, eigvals, x] = fitDynamics(R, epoch_ms, sd_ms, bin_ms, n_pcs)
% fitDynamics  Estimate a linear dynamical system M from neural population data.
%
%   Inputs
%     R        – struct array of trials (same format as mazeJ0918Simple.mat)
%     epoch_ms – [start, stop] ms relative to moveOnsetTime
%     sd_ms    – Gaussian smoothing kernel SD in ms
%     bin_ms   – bin size in ms (assumed 1 ms)
%     n_pcs    – number of PCs to retain
%
%   Outputs
%     M        – n_pcs x n_pcs dynamics matrix
%     eigvals  – eigenvalues of M (column vector)
%     x        – time x n_pcs matrix of PC scores (trial-averaged)

    nTrials  = length(R);
    nUnits   = length(R(1).unit);
    t_start  = epoch_ms(1);
    t_stop   = epoch_ms(2);
    nBins    = t_stop - t_start + 1;   % number of 1-ms bins

    %% ── Part A: Smooth, trial-average, PCA ───────────────────────────────────

    % Build Gaussian smoothing kernel
    kHalf  = ceil(3 * sd_ms);
    kTimes = -kHalf : kHalf;
    kernel = exp(-0.5 * (kTimes / sd_ms).^2);
    kernel = kernel / sum(kernel);

    % Accumulate spike counts in epoch, trial by trial
    % psth: nUnits x nBins x nTrials
    psth = zeros(nUnits, nBins, nTrials);

    for tr = 1:nTrials
        t0 = R(tr).moveOnsetTime;   % alignment event

        for u = 1:nUnits
            spk   = R(tr).unit(u).spikeTimes;
            % Bin into 1-ms bins relative to moveOnsetTime
            edges = (t0 + t_start - 0.5) : bin_ms : (t0 + t_stop + 0.5);
            counts = histcounts(spk, edges);   % 1 x nBins
            % Smooth
            psth(u, :, tr) = conv(counts, kernel, 'same');
        end
    end

    % Trial-average  →  nUnits x nBins
    meanPSTH = mean(psth, 3);   % average over trials

    % PCA on the nBins x nUnits matrix (observations = time points)
    % coeff: nUnits x n_pcs  |  score: nBins x n_pcs
    [~, x, ~] = pca(meanPSTH', 'NumComponents', n_pcs);
    % x is (nBins x n_pcs) — the trajectory in PC space

    %% ── Part B: Estimate x_dot via finite differences ────────────────────────
    % diff along rows (time axis) gives (nBins-1) x n_pcs
    x_dot = diff(x, 1, 1);   % Δx between adjacent time points

    % Use the midpoint states as the x values paired with each x_dot
    x_mid = (x(1:end-1, :) + x(2:end, :)) / 2;   % (nBins-1) x n_pcs

    %% ── Part C: Fit x_dot = M * x  via ordinary linear regression ────────────
    % For each output dimension d:  x_dot(:,d) = x_mid * m_d
    % Solve all dimensions at once:  X_dot = X_mid * M'
    % Normal equations:  M' = (x_mid' * x_mid) \ (x_mid' * x_dot)
    %                    M  = (x_mid \ x_dot)'
    %
    % mldivide gives the least-squares solution column by column.
    M = (x_mid \ x_dot)';   % n_pcs x n_pcs

    %% ── Part D: Eigenvalues ───────────────────────────────────────────────────
    eigvals = eig(M);   % returns column vector of (possibly complex) eigenvalues
end