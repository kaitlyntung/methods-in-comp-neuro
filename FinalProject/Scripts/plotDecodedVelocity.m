function plotDecodedVelocity(decoderModel, R, trialNum)
% Plot actual vs decoded hand velocity over time for a given trial.
% Shows Vx and Vy over time for both ridge and random forest decoders,
% alongside the actual velocity.

binSizeMs = decoderModel.binSizeMs;
nUnits    = length(R(1).unit);

moveOnset = R(trialNum).moveOnsetTime;
moveEnds  = R(trialNum).moveEndsTime;
binEdges  = moveOnset : binSizeMs : moveEnds;
nBins     = length(binEdges) - 1;

% Binned spike times
trialNeural = zeros(nBins, nUnits);
for u = 1:nUnits
    spks = R(trialNum).unit(u).spikeTimes;
    trialNeural(:, u) = histcounts(spks, binEdges) / (binSizeMs / 1000);
end

trialNeural = sqrt(trialNeural);

% Predicted velocities
ridgeVel      = zeros(nBins, 2);
ridgeVel(:,1) = trialNeural * decoderModel.ridge.beta(2:end,1) + decoderModel.ridge.beta(1,1);
ridgeVel(:,2) = trialNeural * decoderModel.ridge.beta(2:end,2) + decoderModel.ridge.beta(1,2);

rfVel = [predict(decoderModel.rf.rfVx, trialNeural), ...
    predict(decoderModel.rf.rfVy, trialNeural)];

% Actual velocity
handTimes  = R(trialNum).HAND.times;
dt         = mean(diff(handTimes));
velX       = gradient(R(trialNum).HAND.X, dt);
velY       = gradient(R(trialNum).HAND.Y, dt);
binCenters = binEdges(1:end-1) + binSizeMs / 2;

actualVel      = zeros(nBins, 2);
actualVel(:,1) = interp1(handTimes, velX, binCenters);
actualVel(:,2) = interp1(handTimes, velY, binCenters);

% Time axis
timeMs = binCenters - moveOnset;

% R^2
r2Ridge = 1 - sum((actualVel(:) - ridgeVel(:)).^2) / sum((actualVel(:) - mean(actualVel(:))).^2);
r2RF    = 1 - sum((actualVel(:) - rfVel(:)).^2)    / sum((actualVel(:) - mean(actualVel(:))).^2);

figure;

% Vx 
subplot(2, 1, 1); hold on;
plot(timeMs, actualVel(:, 1), 'k-', 'LineWidth', 2, 'DisplayName', 'Actual');
plot(timeMs, ridgeVel(:, 1), 'b--', 'LineWidth', 1.5, ...
    'DisplayName', sprintf('Ridge (R^2=%.2f)', r2Ridge));
plot(timeMs, rfVel(:, 1), 'r--', 'LineWidth', 1.5, ...
    'DisplayName', sprintf('RF (R^2=%.2f)', r2RF));
xlabel('Time from movement onset (ms)');
ylabel('V_x (mm/s)');
title('X velocity');
legend('Location', 'best', 'FontSize', 8);
hold off;

% Vy 
subplot(2, 1, 2); hold on;
plot(timeMs, actualVel(:, 2), 'k-', 'LineWidth', 2, 'DisplayName', 'Actual');
plot(timeMs, ridgeVel(:, 2), 'b--', 'LineWidth', 1.5, 'DisplayName', 'Ridge');
plot(timeMs, rfVel(:, 2), 'r--', 'LineWidth', 1.5, 'DisplayName', 'RF');
xlabel('Time from movement onset (ms)');
ylabel('V_y (mm/s)');
title('Y velocity');
legend('Location', 'best', 'FontSize', 8);
hold off;

sgtitle(sprintf('Trial %d — Actual vs Decoded Velocity', trialNum));
end