function plotUMAPByCurve(umapEmbedding, curveDir)
% Plot the UMAP embedding colored by clockwise vs counterclockwise reaches.

isCW       = curveDir == "CW";
isCCW      = curveDir == "CCW";
isStraight = curveDir == "Straight";

figure;
hold on;
scatter(umapEmbedding(isStraight, 1), umapEmbedding(isStraight, 2), ...
    20, [0.7 0.7 0.7], 'filled', 'MarkerFaceAlpha', 0.4, ...
    'DisplayName', 'Straight');
scatter(umapEmbedding(isCCW, 1), umapEmbedding(isCCW, 2), ...
    20, 'b', 'filled', 'MarkerFaceAlpha', 0.6, ...
    'DisplayName', 'Counterclockwise');
scatter(umapEmbedding(isCW, 1), umapEmbedding(isCW, 2), ...
    20, 'r', 'filled', 'MarkerFaceAlpha', 0.6, ...
    'DisplayName', 'Clockwise');

xlabel('Dimension 1');
ylabel('Dimension 2');
title('UMAP — Clockwise vs Counterclockwise Reaches');
legend('Location', 'eastoutside');
hold off;

fprintf('CW: %d, CCW: %d, Straight: %d\n', ...
    sum(isCW), sum(isCCW), sum(isStraight));
end