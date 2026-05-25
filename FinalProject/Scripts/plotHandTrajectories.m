function plotHandTrajectories(R)
% Plot hand trajectories split by maze category
% (straight vs curved) and colored by reach direction.
% Individual trials shown as thin lines, means are the bold lines.
% All trajectories centered at movement onset position. Input is the maze
% data struc array.

% Get Direction Targets
targetXY = reshape([R.targetXY], 2, [])';
angles = atan2d(targetXY(:,2), targetXY(:,1));
angleGroups = round(angles / 45) * 45;
angleGroups(angleGroups == 180) = -180;
angleGroups(angleGroups == 0 & angles < 5)  = -45;   
angleGroups(angleGroups == 0 & angles >= 5) = 45;
uniqueGroups = unique(angleGroups);
nGroups = length(uniqueGroups);

dirLabels = {'Left', 'Lower-Left', 'Lower-Right', 'Upper-Right', 'Up', 'Upper-Left'};
colors = hsv(nGroups);

% Split by maze category
straightMask = [R.mazeCategory] == 'noBarriers';
curvedMask   = [R.mazeCategory] == 'barriers';
categoryMasks = {straightMask, curvedMask};
categoryNames = {'Straight Reaches (No Barriers)', 'Curved Reaches (Barriers)'};

for c = 1:2
    figure;
    hold on;

    meanXAll = cell(nGroups, 1);
    meanYAll = cell(nGroups, 1);

    for g = 1:nGroups
        groupMask = angleGroups == uniqueGroups(g);
        trialIdx  = find(groupMask & categoryMasks{c}');

        if isempty(trialIdx)
            continue
        end

        allX = {};
        allY = {};

        for tr = 1:length(trialIdx)
            trIdx     = trialIdx(tr);
            handX     = R(trIdx).HAND.X;
            handY     = R(trIdx).HAND.Y;
            handTimes = R(trIdx).HAND.times - R(trIdx).moveOnsetTime;

            % Center at movement onset position
            [~, onsetIdx] = min(abs(handTimes));
            handX = handX - handX(onsetIdx);
            handY = handY - handY(onsetIdx);

            % Plot individual trial
            h = plot(handX, handY, 'Color', colors(g,:), ...
                'LineWidth', 0.5, 'HandleVisibility', 'off');
            h.Color(4) = 0.05;

            allX{end+1} = handX;
            allY{end+1} = handY;
        end

        % Calculate mean
        nPoints = 100;
        meanX   = zeros(length(allX), nPoints);
        meanY   = zeros(length(allY), nPoints);

        for tr = 1:length(allX)
            tNorm        = linspace(0, 1, length(allX{tr}));
            tInterp      = linspace(0, 1, nPoints);
            meanX(tr, :) = interp1(tNorm, allX{tr}, tInterp);
            meanY(tr, :) = interp1(tNorm, allY{tr}, tInterp);
        end

        meanXAll{g} = mean(meanX, 1);
        meanYAll{g} = mean(meanY, 1);
    end

    % Plot means
    for g = 1:nGroups
        if ~isempty(meanXAll{g})
            plot(meanXAll{g}, meanYAll{g}, 'k-', 'LineWidth', 7, ...
                'HandleVisibility', 'off');
            plot(meanXAll{g}, meanYAll{g}, 'Color', colors(g,:), ...
                'LineWidth', 5, 'DisplayName', dirLabels{g});
        end
    end

    % Movement onset marker
    plot(0, 0, 'k+', 'MarkerSize', 12, 'LineWidth', 2, 'HandleVisibility', 'off');

    xlabel('X position');
    ylabel('Y position');
    title(categoryNames{c});
    legend('Location', 'eastoutside', 'FontSize', 8);
    axis equal;
    hold off;
end
end