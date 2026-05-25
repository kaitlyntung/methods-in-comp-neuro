function plotDecodedTrajectory(decoderModel, R, trialNum)
% Plot actual vs decoded hand trajectory for a chosen trial.
% Uses trained ridge and random forest models from trainDecoder.
% Trajectories are normalized to unit path length.

binSizeMs = decoderModel.binSizeMs;
nUnits    = length(R(1).unit);

moveOnset = R(trialNum).moveOnsetTime;
moveEnds  = R(trialNum).moveEndsTime;
binEdges  = moveOnset : binSizeMs : moveEnds;
nBins     = length(binEdges) - 1;

trialNeural = zeros(nBins, nUnits);
for u = 1:nUnits
    spks = R(trialNum).unit(u).spikeTimes;
    trialNeural(:, u) = histcounts(spks, binEdges) / (binSizeMs / 1000);
end

trialNeural = sqrt(trialNeural);

ridgeVel      = zeros(nBins, 2);
ridgeVel(:,1) = trialNeural * decoderModel.ridge.beta(2:end,1) + decoderModel.ridge.beta(1,1);
ridgeVel(:,2) = trialNeural * decoderModel.ridge.beta(2:end,2) + decoderModel.ridge.beta(1,2);

rfVel = [predict(decoderModel.rf.rfVx, trialNeural), ...
    predict(decoderModel.rf.rfVy, trialNeural)];

handTimes  = R(trialNum).HAND.times;
dt         = mean(diff(handTimes));
velX       = gradient(R(trialNum).HAND.X, dt);
velY       = gradient(R(trialNum).HAND.Y, dt);
binCenters = binEdges(1:end-1) + binSizeMs / 2;

actualVel      = zeros(nBins, 2);
actualVel(:,1) = interp1(handTimes, velX, binCenters);
actualVel(:,2) = interp1(handTimes, velY, binCenters);

% Trial R^2
r2Ridge = 1 - sum((actualVel(:) - ridgeVel(:)).^2) / sum((actualVel(:) - mean(actualVel(:))).^2);
r2RF    = 1 - sum((actualVel(:) - rfVel(:)).^2)    / sum((actualVel(:) - mean(actualVel(:))).^2);

dtSec      = binSizeMs / 1000;
actualTraj = normalizeTrajectory(cumsum(actualVel * dtSec, 1));
ridgeTraj  = normalizeTrajectory(cumsum(ridgeVel  * dtSec, 1));
rfTraj     = normalizeTrajectory(cumsum(rfVel     * dtSec, 1));

figure;
subplot(1,2,1); hold on;
plot(actualTraj(:,1), actualTraj(:,2), 'k-',  'LineWidth', 2,   'DisplayName', 'Actual');
plot(ridgeTraj(:,1),  ridgeTraj(:,2),  'b--', 'LineWidth', 1.5, 'DisplayName', sprintf('Ridge R²=%.2f', r2Ridge));
title('Ridge Regression'); legend; axis equal; hold off;

subplot(1,2,2); hold on;
plot(actualTraj(:,1), actualTraj(:,2), 'k-',  'LineWidth', 2,   'DisplayName', 'Actual');
plot(rfTraj(:,1),     rfTraj(:,2),     'r--', 'LineWidth', 1.5, 'DisplayName', sprintf('RF R²=%.2f', r2RF));
xlabel('X position'); ylabel('Y position');
title('Random Forest'); legend; axis equal; hold off;

sgtitle(sprintf('Trial %d — Reconstructed Hand Trajectory', trialNum));
end