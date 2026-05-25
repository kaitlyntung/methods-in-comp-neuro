function plotNLDimReduction(pcaData, conditionLabels)
% Apply t-SNE and UMAP to PCA-reduced neural data and plot embeddings
% colored by reach condition. Compute pairwise distance matrix
% correlation to compare the two methods.
%
% Inputs:
%   pcaData          - nTrials x nDims PCA scores
%   conditionLabels  - nTrials x 1 vector of direction group labels

uniqueConditions = unique(conditionLabels);
nConditions = length(uniqueConditions);
colors = hsv(nConditions);

dirLabels = {'Left', 'Lower-Left', 'Lower-Right', 'Upper-Right', 'Up', 'Upper-Left'};

% Run t-SNE with PCA initialization
rng(1);
tSNEEmbedding = tsne(pcaData, 'InitialY', pcaData(:,1:2), ...
    'Perplexity', 30);

% Run UMAP
rng(1);
umapEmbedding = umap(pcaData, Reproducible=true);

% Plot side by side
figure;
embeddings = {tSNEEmbedding, umapEmbedding};
embeddingNames = {'t-SNE', 'UMAP'};

for e = 1:2
    subplot(1, 2, e);
    hold on;

    for c = 1:nConditions
        condMask = conditionLabels == uniqueConditions(c);
        scatter(embeddings{e}(condMask, 1), embeddings{e}(condMask, 2), ...
            20, colors(c,:), 'filled', 'MarkerFaceAlpha', 0.6, ...
            'DisplayName', dirLabels{c});
    end

    xlabel('Dimension 1');
    ylabel('Dimension 2');
    title(embeddingNames{e});
    hold off;
end

legend('Location', 'eastoutside', 'FontSize', 7);
sgtitle('M1/PMd Population Activity — Nonlinear Dimensionality Reduction');

% Pairwise distance matrices:
tSNEDists = pdist(tSNEEmbedding);
umapDists = pdist(umapEmbedding);
r = corr(tSNEDists', umapDists');
fprintf('Pairwise distance matrix correlation (t-SNE vs UMAP): r = %.3f\n', r);
end