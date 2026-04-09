%% 2.1 Visualization: scatter plots
% Part A
colors = {'b', 'r', 'k'};
titles = {'Type 1', 'Type 2', 'Type 3'};

for neuron = 1:3
    figure;
    location = data(neuron).locsXY;
    
    plot(location(:,1), location(:,2), 'o', ...
        'MarkerEdgeColor', colors{neuron}, ...
        'MarkerFaceColor', colors{neuron}, ...
        'MarkerSize', 4);
    hold on;
    
    fov = data(neuron).fovSizeXY;
    rectangle('Position', [0, 0, fov(1), fov(2)], ...
              'EdgeColor', 'k', 'LineWidth', 2);
    
    axis equal;
    xlim([0 fov(1)]);
    ylim([0 fov(2)]);
    
    xlabel('X position');
    ylabel('Y position');
    title([titles{neuron} ' Scatter Plot']);
end

%% 2.2 Nearest-Neighbor Plots
% Part A
function nnDist = nearestNeighborDistances(locsXY)
    D = pdist2(locsXY, locsXY);
    D(1:size(D,1)+1:end) = inf;
    nnDist = min(D, [], 2);
end

nnDists = cell(1,3);
medians = zeros(1,3);

for i = 1:3
    locs = data(i).locsXY;
    
    nnDists{i} = nearestNeighborDistances(locs);
    medians(i) = median(nnDists{i});
end
allDists = vertcat(nnDists{:});
edges = linspace(0, max(allDists), 30);
colors = {'b','r','k'};

for i = 1:3
    figure;
    
    histogram(nnDists{i}, edges, 'FaceColor', colors{i});
    hold on;
    plot(medians(i), 0, 'rx', 'MarkerSize', 12, 'LineWidth', 2);
    xlabel('Nearest Neighbor Distance');
    ylabel('Count');
    title(['Neuron Type ' num2str(i) ' Histogram']);
    xlim([edges(1) edges(end)]);
end

for i = 1:3
    locs = data(i).locsXY;
    D = pdist2(locs, locs);
    D(1:size(D,1)+1:end) = inf;
    [~, nnIdx] = min(D, [], 2);

    figure;
    hold on;
    plot(locs(:,1), locs(:,2), 'o', ...
        'MarkerFaceColor', colors{i}, ...
        'MarkerEdgeColor', colors{i});
    for j = 1:size(locs,1)
        x = [locs(j,1), locs(nnIdx(j),1)];
        y = [locs(j,2), locs(nnIdx(j),2)];
        
        plot(x, y, colors{i});
    end
    
    fov = data(i).fovSizeXY;
    rectangle('Position', [0 0 fov(1) fov(2)], 'EdgeColor', 'k');
    axis equal;
    title(['Neuron Type ' num2str(i) ' with Nearest Neighbors']);
end

%% Part 2.3: Testing for Clustery-ness
% Part A
function locs = generateRandomLocs(n, fov)
    x = rand(n,1) * fov(1);
    y = rand(n,1) * fov(2);
    locs = [x y];
end

% Part B
numRuns = 1000;
mcMedians = zeros(numRuns,1);
locs_real = data(1).locsXY;
fov = data(1).fovSizeXY;
n = size(locs_real,1);

real_nn = nearestNeighborDistances(locs_real);
real_median = median(real_nn);

for i = 1:numRuns
    % Generate random neurons
    locs_rand = generateRandomLocs(n, fov);
    % Compute NN distances
    nn_rand = nearestNeighborDistances(locs_rand);
    % Store median
    mcMedians(i) = median(nn_rand);
end

% Part C
p_lower = mean(mcMedians <= real_median);
p_upper = mean(mcMedians >= real_median);
% choose the smaller tail
p_one_tailed = min(p_lower, p_upper);
% convert to two-tailed
p_two_tailed = 2 * p_one_tailed;
if real_median < mean(mcMedians)
    direction = 'more clustered';
else
    direction = 'less clustered';
end

% Part D
fprintf('Cell type 1 is %s than random, p = %.3f\n', ...
        direction, p_two_tailed);

%% Repeat for all three neuron types
numRuns = 1000;

for neuron = 1:3
    
    locs_real = data(neuron).locsXY;
    fov = data(neuron).fovSizeXY;
    n = size(locs_real,1);
    
    real_nn = nearestNeighborDistances(locs_real);
    real_median = median(real_nn);
    
    mcMedians = zeros(numRuns,1);
    
    for i = 1:numRuns
        locs_rand = generateRandomLocs(n, fov);
        nn_rand = nearestNeighborDistances(locs_rand);
        mcMedians(i) = median(nn_rand);
    end
    
    p_lower = mean(mcMedians <= real_median);
    p_upper = mean(mcMedians >= real_median);
    
    p_two = 2 * min(p_lower, p_upper);
    
    if real_median < mean(mcMedians)
        direction = 'more clustered';
    else
        direction = 'less clustered';
    end
    
    fprintf('Cell type %d is %s than random, p = %.3f\n', ...
        neuron, direction, p_two);
end