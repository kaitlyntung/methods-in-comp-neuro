function plotNeurometricFunctions(R, threshold)
% Plot neurometric functions for all of the motor units using ROC-based Ideal
% Observer Analysis during the pre-movement period, split by brain area.
% Each neuron's discriminability (AUROC) is plotted as a function of
% angular distance from its preferred reach direction. Inputs: maze data
% and pre movement window length(ms)/threshold.

nUnits = length(R(1).unit);
whichArray = R(1).whichArray;

% Get Reach Targets
targetXY = reshape([R.targetXY], 2, [])';
angles = atan2d(targetXY(:,2), targetXY(:,1));
angleGroups = round(angles / 45) * 45;
angleGroups(angleGroups == 180) = -180;
angleGroups(angleGroups == 0 & angles < 5)  = -45;
angleGroups(angleGroups == 0 & angles >= 5) = 45;
uniqueGroups = unique(angleGroups);
nGroups = length(uniqueGroups);

% Get differences
allAngularDiffs = zeros(nGroups, nGroups);
for g1 = 1:nGroups
    for g2 = 1:nGroups
        diff = abs(uniqueGroups(g1) - uniqueGroups(g2));
        if diff > 180
            diff = 360 - diff;
        end
        allAngularDiffs(g1, g2) = diff;
    end
end
uniqueDiffs = unique(allAngularDiffs(:));
nDiffs = length(uniqueDiffs);

% Pre-allocate matrix
aurocByDiff = zeros(nUnits, nDiffs);

for u = 1:nUnits

    meanDelayFR = zeros(nGroups, 1);
    spkCountPerGroup = cell(nGroups, 1);

    for g = 1:nGroups
        groupTrials = find(angleGroups == uniqueGroups(g));
        frs = [];

        for tr = 1:length(groupTrials)
            trIdx     = groupTrials(tr);
            spks      = R(trIdx).unit(u).spikeTimes;
            moveOnset = R(trIdx).moveOnsetTime;

            % Filtering
            trialAvailable = moveOnset - R(trIdx).targetAppearsTime;
            if trialAvailable > threshold
                alignedSpks = spks - moveOnset;
                preMovSpks  = sum(alignedSpks >= -threshold & alignedSpks < 0);
                frs(end + 1) = preMovSpks / (threshold / 1000);
            end
        end

        meanDelayFR(g)        = mean(frs);
        spkCountPerGroup{g}   = frs;
    end

    % Find preferred direction
    [~, prefIdx] = max(meanDelayFR);
    prefSpkCount = spkCountPerGroup{prefIdx}';

    % Compute AUROC for each direction vs preferred
    aurocPerGroup  = zeros(nGroups, 1);
    angDiffPerGroup = zeros(nGroups, 1);

    for g = 1:nGroups
        compSpkCount = spkCountPerGroup{g}';

        if g == prefIdx
            aurocPerGroup(g) = 0.5;
        else
            labels = [ones(length(prefSpkCount), 1); zeros(length(compSpkCount), 1)];
            scores = [prefSpkCount; compSpkCount];
            [~, ~, ~, auroc] = perfcurve(labels, scores, 1);
            aurocPerGroup(g) = auroc;
        end

        % Angular difference
        angDiff = abs(uniqueGroups(g) - uniqueGroups(prefIdx));
        if angDiff > 180
            angDiff = 360 - angDiff;
        end
        angDiffPerGroup(g) = angDiff;
    end

    % Average AUROC 
    for d = 1:nDiffs
        mask = angDiffPerGroup == uniqueDiffs(d);
        aurocByDiff(u, d) = mean(aurocPerGroup(mask));
    end
end

% Split by area
pmdMask   = whichArray == 1;
m1Mask    = whichArray == 2;
areaMasks = {pmdMask, m1Mask};
areaNames = {'PMd', 'M1'};
meanColors = {'b', 'r'};

figure;
for a = 1:2
    subplot(1, 2, a);
    hold on;

    % Individual neurons
    areaUnits = find(areaMasks{a});
    for i = 1:length(areaUnits)
        u = areaUnits(i);
        validMask = ~isnan(aurocByDiff(u, :));
        plot(uniqueDiffs(validMask), aurocByDiff(u, validMask), ...
            'Color', [0.75 0.75 0.75 0.15], 'LineWidth', 0.5, ...
            'HandleVisibility', 'off');
    end

    % Population mean
    popMean = nanmean(aurocByDiff(areaMasks{a}, :), 1);
    plot(uniqueDiffs, popMean, meanColors{a}, 'LineWidth', 3, ...
        'DisplayName', sprintf('%s mean', areaNames{a}));

    % Threshold line
    yline(0.75, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');

    xlim([0 180]);
    xticks(0:45:180);
    ylim([0.4 1]);
    xlabel('\Delta\theta (deg)');
    ylabel('AUROC');
    title(sprintf('%s (n=%d)', areaNames{a}, sum(areaMasks{a})));
    legend('Location', 'northwest');
    hold off;
end

sgtitle('Neurometric Functions — Pre-Movement Activity');

% How many neurons crossed threshold?
threshold75 = 0.75;

pmdUnits = find(R(1).whichArray == 1);
m1Units  = find(R(1).whichArray == 2);

pmdAuroc = aurocByDiff(pmdUnits, :);
m1Auroc  = aurocByDiff(m1Units,  :);

pmdCross = sum(any(pmdAuroc > threshold75, 2));
m1Cross  = sum(any(m1Auroc  > threshold75, 2));

fprintf('PMd: %d/%d neurons (%.0f%%) exceed AUROC 0.75 at any distance\n', ...
    pmdCross, length(pmdUnits), 100*pmdCross/length(pmdUnits));
fprintf('M1:  %d/%d neurons (%.0f%%) exceed AUROC 0.75 at any distance\n', ...
    m1Cross,  length(m1Units),  100*m1Cross/length(m1Units));
end