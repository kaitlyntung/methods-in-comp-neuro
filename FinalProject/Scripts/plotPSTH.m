function plotPSTH(R, unitNum, preMovMs, postMovMs, binSizeMs)
% Plot a PSTH for a motor unit aligned to movement onset,
% averaged by reach direction group with SEM error bounds. Inputs include:
% maze data set, unit to plot, ms before movement onset, ms after movement
% onset, and bin size.

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

% Time axis
binEdges = -preMovMs : binSizeMs : postMovMs;
binCenters = binEdges(1:end-1) + binSizeMs/2;
nBins = length(binCenters);

% Brain area for this unit
arrayLabel = R(1).arrays{R(1).whichArray(unitNum)};

figure;
hold on;

for g = 1:nGroups
    % Direction group trials
    groupMask = angleGroups == uniqueGroups(g);
    groupTrials = find(groupMask);
    nTrials = length(groupTrials);

    % Preallocate Matrix
    spikeCounts = zeros(nTrials, nBins);

    for tr = 1:nTrials
        trIdx = groupTrials(tr);

        % Align spike times to movement onset
        spks = R(trIdx).unit(unitNum).spikeTimes;
        moveOnset = R(trIdx).moveOnsetTime;
        alignedSpks = spks - moveOnset;

        % Bin spikes
        spikeCounts(tr, :) = histcounts(alignedSpks, binEdges);
    end

    % Convert to firing rate in spikes/s
    firingRate = spikeCounts / (binSizeMs / 1000);

    % Mean and SEM across trials
    meanFR = mean(firingRate, 1);
    semFR = std(firingRate, 0, 1) / sqrt(nTrials);
    meanFR = smoothdata(meanFR, 'gaussian', 7);
    semFR  = smoothdata(semFR,  'gaussian', 7);

    % SEM error bounds
    upper = meanFR + semFR;
    lower = meanFR - semFR;
    fill([binCenters, fliplr(binCenters)], ...
        [upper, fliplr(lower)], ...
        colors(g,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none', ...
        'HandleVisibility', 'off');

    % Mean line
    plot(binCenters, meanFR, 'Color', colors(g,:), 'LineWidth', 1.5, ...
        'DisplayName', dirLabels{g});
end

% Movement onset line
xline(0, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');

xlabel('Time relative to movement onset (ms)');
ylabel('Firing Rate (spikes/s)');
xlim([-preMovMs, postMovMs]);
title(sprintf('Unit %d (%s) PSTH by Reach Direction', unitNum, arrayLabel));
legend('Location', 'eastoutside', 'FontSize', 8);
hold off;
end