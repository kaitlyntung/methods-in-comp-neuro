function ridgeModel = trainRidgeDecoder(neuralData, velocityData, nFolds)
% Train a ridge regression decoder to predict hand velocity (Vx, Vy)
% from firing rates. Finds optimal lambda via k-fold CV,
% then trains final model on all data.

lambdaValues = logspace(-4, 4, 50);
nLambdas     = length(lambdaValues);
nTrials      = length(neuralData);

cv        = cvpartition(nTrials, 'KFold', nFolds);
lambdaMSE = zeros(nLambdas, 1);

for fold = 1:nFolds
    neuralTrain = vertcat(neuralData{training(cv, fold)});
    velTrain    = vertcat(velocityData{training(cv, fold)});
    neuralTest  = vertcat(neuralData{test(cv, fold)});
    velTest     = vertcat(velocityData{test(cv, fold)});

    lambdaWeightsVx = ridge(velTrain(:,1), neuralTrain, lambdaValues, 0);
    lambdaWeightsVy = ridge(velTrain(:,2), neuralTrain, lambdaValues, 0);

    for lIdx = 1:nLambdas
        predVx = neuralTest * lambdaWeightsVx(2:end, lIdx) + lambdaWeightsVx(1, lIdx);
        predVy = neuralTest * lambdaWeightsVy(2:end, lIdx) + lambdaWeightsVy(1, lIdx);
        lambdaMSE(lIdx) = lambdaMSE(lIdx) + ...
            (mean((velTest(:,1) - predVx).^2) + mean((velTest(:,2) - predVy).^2)) / 2;
    end
end

lambdaMSE    = lambdaMSE / nFolds;
[~, bestIdx] = min(lambdaMSE);
bestLambda   = lambdaValues(bestIdx);

% CV and R^2
foldR2 = zeros(nFolds, 2);
for fold = 1:nFolds
    neuralTrain = vertcat(neuralData{training(cv, fold)});
    velTrain    = vertcat(velocityData{training(cv, fold)});
    neuralTest  = vertcat(neuralData{test(cv, fold)});
    velTest     = vertcat(velocityData{test(cv, fold)});

    wVx = ridge(velTrain(:,1), neuralTrain, bestLambda, 0);
    wVy = ridge(velTrain(:,2), neuralTrain, bestLambda, 0);

    predVx = neuralTest * wVx(2:end) + wVx(1);
    predVy = neuralTest * wVy(2:end) + wVy(1);

    foldR2(fold,1) = 1 - sum((velTest(:,1) - predVx).^2) / sum((velTest(:,1) - mean(velTest(:,1))).^2);
    foldR2(fold,2) = 1 - sum((velTest(:,2) - predVy).^2) / sum((velTest(:,2) - mean(velTest(:,2))).^2);
end

r2 = mean(foldR2, 1);
fprintf('Ridge: best lambda: %.4f , CV R²: Vx=%.3f, Vy=%.3f\n', bestLambda, r2(1), r2(2));

% Train model on all data
neuralAll = vertcat(neuralData{:});
velAll    = vertcat(velocityData{:});

finalWeightsVx = ridge(velAll(:,1), neuralAll, bestLambda, 0);
finalWeightsVy = ridge(velAll(:,2), neuralAll, bestLambda, 0);

ridgeModel.beta         = [finalWeightsVx, finalWeightsVy];
ridgeModel.lambda       = bestLambda;
ridgeModel.r2           = r2;
ridgeModel.lambdaMSE    = lambdaMSE;
ridgeModel.lambdaValues = lambdaValues;
end