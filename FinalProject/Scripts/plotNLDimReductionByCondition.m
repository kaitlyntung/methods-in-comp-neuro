function [tSNEEmbedding, umapEmbedding] = plotNLDimReductionByCondition(pcaData, conditionLabels, isStraight)
% Apply t-SNE and UMAP to PCA-reduced neural data and plot embeddings
% colored by reach direction, with maze condition encoded by marker outline
% 
% Inputs:
% pcaData          - nTrials x nDims PCA scores
% conditionLabels  - nTrials x 1 vector of direction group labels
% isStraight = nTrials x 1 vector logical of straight vs. not

uniqueConditions = unique(conditionLabels);
nConditions = length(uniqueConditions);
colors = hsv(nConditions);

dirLabels = {'Left', 'Lower-Left', 'Lower-Right', 'Upper-Right', 'Up', 'Upper-Left'};

% t-SNE and UMAP
rng(1);
tSNEEmbedding = tsne(pcaData, 'InitialY', pcaData(:,1:2), ...
    'Perplexity', 30);
rng(1);
umapEmbedding = umap(pcaData, Reproducible=true);

figure;
embeddings = {tSNEEmbedding, umapEmbedding};
embeddingNames = {'t-SNE', 'UMAP'};

for e = 1:2
    subplot(1, 2, e);
    hold on;

    % Direction
    for c = 1:nConditions
        dirMask = conditionLabels == uniqueConditions(c);

        % Straight
        straightMask = dirMask & isStraight;
        scatter(embeddings{e}(straightMask, 1), embeddings{e}(straightMask, 2), ...
            25, colors(c,:), 'filled', 'MarkerFaceAlpha', 0.6, ...
            'DisplayName', dirLabels{c});

        % Curved
        curvedMask = dirMask & ~isStraight;
        scatter(embeddings{e}(curvedMask, 1), embeddings{e}(curvedMask, 2), ...
            25, colors(c,:), 'filled', 'MarkerEdgeColor', 'k', ...
            'LineWidth', 0.6, 'MarkerFaceAlpha', 0.6, ...
            'HandleVisibility', 'off');
    end

    % Dummy entries for legend
    scatter(NaN, NaN, 40, [0.5 0.5 0.5], 'filled', 'MarkerFaceAlpha', 0.6, ...
        'DisplayName', 'Straight');
    scatter(NaN, NaN, 40, [0.5 0.5 0.5], 'filled', 'MarkerEdgeColor', 'k', ...
        'LineWidth', 0.6, 'MarkerFaceAlpha', 0.6, 'DisplayName', 'Curved');

    xlabel('Dimension 1');
    ylabel('Dimension 2');
    title(embeddingNames{e});
    if e == 2
        legend('Location', 'eastoutside', 'FontSize', 7);
    end
    hold off;
end

sgtitle('M1/PMd Population Activity — by Direction and Maze Condition');

% Pairwise distance correlation
tSNEDists = pdist(tSNEEmbedding);
umapDists = pdist(umapEmbedding);
r = corr(tSNEDists', umapDists');
fprintf('Pairwise distance matrix correlation (t-SNE vs UMAP, by-condition): r = %.3f\n', r);
fprintf('Straight trials: %d, Curved trials: %d\n', ...
    sum(isStraight), sum(~isStraight));
end