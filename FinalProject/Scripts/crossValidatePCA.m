function [optimalDims, pcaData] = crossValidatePCA(firingRateMatrix, maxDims)
% Double Cross Validation to find the optimal dimensionality.
% Both observations (trials) and variables (neurons) are split into
% k folds. For each pair, we fit PCA on the
% training quadrant and evaluate the reconstruction for the held out test quadrant.

rng(1);
nTrials = size(firingRateMatrix, 1);
nUnits  = size(firingRateMatrix, 2);
nFolds  = 10;

% Trial folds
trialCv = cvpartition(nTrials, 'KFold', nFolds);

% Neuron fold
neuronPerm  = randperm(nUnits);
neuronFolds = mod(0:nUnits-1, nFolds) + 1;
neuronFoldId = zeros(1, nUnits);
neuronFoldId(neuronPerm) = neuronFolds;

% Pre-allocate error matrix
mse = zeros(nFolds, nFolds, maxDims);

for trialFold = 1:nFolds
    trainTrialMask = training(trialCv, trialFold);
    testTrialMask  = test(trialCv, trialFold);

    for neuronFold = 1:nFolds
        yNeurons = find(neuronFoldId == neuronFold);
        xNeurons = find(neuronFoldId ~= neuronFold);

        % Four quadrants
        xTrain = firingRateMatrix(trainTrialMask, xNeurons);
        yTrain = firingRateMatrix(trainTrialMask, yNeurons);
        xTest  = firingRateMatrix(testTrialMask,  xNeurons);
        yTest  = firingRateMatrix(testTrialMask,  yNeurons);

        % Subtract training means
        xTrainMean = mean(xTrain, 1);
        yTrainMean = mean(yTrain, 1);
        xTrainC = xTrain - xTrainMean;
        xTestC  = xTest  - xTrainMean;
        yTrainC = yTrain - yTrainMean;

        % PCA on xTrainC, request all maxDims components at once
        [V, ~] = pca(xTrainC, 'NumComponents', maxDims);

        % Project xTrain and xTest once
        lowDTrain = xTrainC * V;
        lowDTest  = xTestC  * V;

        for d = 1:maxDims
            % Regress yTrain (centered) on top-d PC scores of xTrain
            scoresTrain = lowDTrain(:, 1:d);
            beta = (scoresTrain' * scoresTrain) \ (scoresTrain' * yTrainC);

            % Predict yTest from top-d PC scores of xTest, add back yTrainMean
            scoresTest = lowDTest(:, 1:d);
            yPred = scoresTest * beta + yTrainMean;

            mse(trialFold, neuronFold, d) = mean((yTest - yPred).^2, 'all');
        end
    end
end

% Average across both fold dimensions
meanError = squeeze(mean(mse, [1 2]));

% Plot error curve
figure;
plot(1:maxDims, meanError, 'k', 'LineWidth', 2);
xlabel('Number of PCA Components');
ylabel('Mean Squared Error');
title('Double Cross-Validated PCA — Optimal Dimensionality');

% Optimal D: minimum of the mean error curve
[~, optimalDims] = min(meanError);
xline(optimalDims, 'r--', 'LineWidth', 1.5, ...
    'Label', sprintf('Optimal: %d dims', optimalDims));

% Final PCA on all data with optimal dims
[~, pcaData] = pca(firingRateMatrix, 'NumComponents', optimalDims);

fprintf('Optimal dimensionality: %d\n', optimalDims);
end