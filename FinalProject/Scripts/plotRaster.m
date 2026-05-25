function plotRaster(R, unitNum, preMovMs, postMovMs)
% Plot a raster for a single motor unit aligned to movement onset,
% with trials sorted and colored by reach direction group. Inputs include:
% maze data set, unit number, pre and post movement windows.

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

% Brain area
arrayLabel = R(1).arrays{R(1).whichArray(unitNum)};

figure;
hold on;

trialOffset = 0;

for g = 1:nGroups
    groupMask = angleGroups == uniqueGroups(g);
    groupTrials = find(groupMask);

    % Dummy line for legend
    plot(NaN, NaN, 'Color', colors(g,:), 'LineWidth', 2, ...
        'DisplayName', dirLabels{g});

    for tr = 1:length(groupTrials)
        trIdx = groupTrials(tr);

        % Align spike times to movement onset
        spks = R(trIdx).unit(unitNum).spikeTimes;
        moveOnset = R(trIdx).moveOnsetTime;
        alignedSpks = spks - moveOnset;

        % Window
        windowMask = alignedSpks >= -preMovMs & alignedSpks <= postMovMs;
        alignedSpks = alignedSpks(windowMask);

        if isempty(alignedSpks)
            continue
        end

        % Plot spike ticks
        y1 = trialOffset + tr;
        y2 = y1 + 1;
        spksRep = repmat(alignedSpks', [2 1]);
        ypos = repmat([y1; y2], [1 length(alignedSpks)]);
        plot(spksRep, ypos, 'Color', colors(g,:), ...
            'LineWidth', 1.5, 'HandleVisibility', 'off');
    end

    trialOffset = trialOffset + length(groupTrials);
end

% Movement onset line
xline(0, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');

set(gca, 'YDir', 'reverse', 'YLim', [0 trialOffset]);
xlabel('Time relative to movement onset (ms)');
ylabel('Trial');
xlim([-preMovMs, postMovMs]);
title(sprintf('Unit %d (%s) Raster Plot', unitNum, arrayLabel));
legend('Location', 'eastoutside', 'FontSize', 8);
hold off;
end