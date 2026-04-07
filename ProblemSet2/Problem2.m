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
    
    xmin = min(location(:,1));
    xmax = max(location(:,1));
    ymin = min(location(:,2));
    ymax = max(location(:,2));
    rectangle('Position', [xmin, ymin, xmax-xmin, ymax-ymin], ...
              'EdgeColor', 'k', 'LineWidth', 2);
    
    axis equal;
    xlabel('X position');
    ylabel('Y position');
    title([titles{neuron} ' Scatter Plot']);
    
end

%% 2.2 Nearest-Neighbor Plots
% Part A
