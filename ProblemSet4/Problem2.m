%% Part A
% predictors X1, outcomes y1

% part i
cv = cvpartition(length(y1), 'KFold', 5);
num_folds = cv.NumTestSets;
error = zeros(num_folds, 1);
VAF_train = zeros(num_folds, 1);
VAF_test  = zeros(num_folds, 1);

figure;
for fold = 1:num_folds
    trainInd = training(cv, fold);
    testInd = test(cv, fold);

    Xtrain = X1(trainInd, :);
    ytrain = y1(trainInd);
    Xtest = X1(testInd, :);
    ytest = y1(testInd);

    Xtrain = [ones(size(Xtrain, 1), 1), Xtrain];
    Xtest = [ones(size(Xtest, 1), 1), Xtest];

    b = Xtrain \ ytrain;
    ypred = Xtest * b;

    % part iii
    % Compute VAF (training)
    VAF_train(fold) = 1 - sum((ytrain - ypred_train).^2) / sum(ytrain.^2);
    
    % Compute VAF (test)
    VAF_test(fold) = 1 - sum((ytest - ypred_test).^2) / sum(ytest.^2);
    
    % Print per fold
    fprintf('Fold %d: Train VAF = %.3f, Test VAF = %.3f\n', ...
            fold, VAF_train(fold), VAF_test(fold));

    % part ii
    subplot(2,3,fold);
    scatter(ytest, ypred, 'filled'); hold on;
    
    minVal = min([ytest; ypred]);
    maxVal = max([ytest; ypred]);
    plot([minVal maxVal], [minVal maxVal], 'k-', 'LineWidth', 1.5);
    
    xlabel('True y');
    ylabel('Predicted y');
    title(['Fold ' num2str(fold)]);
    
    error(fold) = mean((ytest - ypred).^2);
end

mean_error = mean(error);
disp(mean_error);
