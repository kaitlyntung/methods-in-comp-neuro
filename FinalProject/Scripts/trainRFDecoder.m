function rfModel = trainRFDecoder(neuralData, velocityData, nFolds)
% Train a random forest decoder to predict hand velocity (Vx, Vy)
% from firing rates. Gets R^2 from k-fold CV, and 
% then trains final models on all data.

nTrials = length(neuralData);
nTrees  = 100;
cv      = cvpartition(nTrials, 'KFold', nFolds);
foldR2  = zeros(nFolds, 2);

for fold = 1:nFolds
    neuralTrain = vertcat(neuralData{training(cv, fold)});
    velTrain    = vertcat(velocityData{training(cv, fold)});
    neuralTest  = vertcat(neuralData{test(cv, fold)});
    velTest     = vertcat(velocityData{test(cv, fold)});

    rfVx  = TreeBagger(nTrees, neuralTrain, velTrain(:,1), 'Method', 'regression', 'OOBPrediction', 'off');
    rfVy  = TreeBagger(nTrees, neuralTrain, velTrain(:,2), 'Method', 'regression', 'OOBPrediction', 'off');
    yPred = [predict(rfVx, neuralTest), predict(rfVy, neuralTest)];

    foldR2(fold,1) = 1 - sum((velTest(:,1) - yPred(:,1)).^2) / sum((velTest(:,1) - mean(velTest(:,1))).^2);
    foldR2(fold,2) = 1 - sum((velTest(:,2) - yPred(:,2)).^2) / sum((velTest(:,2) - mean(velTest(:,2))).^2);
end

meanR2 = mean(foldR2, 1);
fprintf('Random Forest — R²: Vx=%.3f, Vy=%.3f\n', meanR2(1), meanR2(2));

% Train model on all data
neuralAll = vertcat(neuralData{:});
velAll    = vertcat(velocityData{:});

finalRFVx = TreeBagger(nTrees, neuralAll, velAll(:,1), 'Method', 'regression', 'OOBPrediction', 'off');
finalRFVy = TreeBagger(nTrees, neuralAll, velAll(:,2), 'Method', 'regression', 'OOBPrediction', 'off');

rfModel.rfVx   = finalRFVx;
rfModel.rfVy   = finalRFVy;
rfModel.nTrees = nTrees;
rfModel.r2     = meanR2;
end