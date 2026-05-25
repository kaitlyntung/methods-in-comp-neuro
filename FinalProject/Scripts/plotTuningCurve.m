function plotTuningCurve(R)
% Plot population tuning curves for all motor units during the movement
% period, normalized to each neuron's peak firing rate. The individual neurons are plotted as transparent 
% lines, and the population as the bold line. Input: maze data set.

nUnits = length(R(1).unit);

% Get Direction Target
targetXY = reshape([R.targetXY], 2, [])';
angles = atan2d(targetXY(:,2), targetXY(:,1));
angleGroups = round(angles / 45) * 45;
angleGroups(angleGroups == 180) = -180;
angleGroups(angleGroups == 0 & angles < 5)  = -45;
angleGroups(angleGroups == 0 & angles >= 5) = 45;
uniqueGroups = unique(angleGroups);
nGroups = length(uniqueGroups);

% Convert to radians
dirRad = deg2rad(uniqueGroups);
dirRadClosed = [dirRad; dirRad(1)];

% Compute FR
allFR = zeros(nUnits, nGroups);

for u = 1:nUnits
    for g = 1:nGroups
        groupTrials = find(angleGroups == uniqueGroups(g));
        frs = zeros(length(groupTrials), 1);

        for tr = 1:length(groupTrials)
            trIdx     = groupTrials(tr);
            spks      = R(trIdx).unit(u).spikeTimes;
            moveOnset = R(trIdx).moveOnsetTime;
            moveEnds  = R(trIdx).moveEndsTime;
            movDur    = moveEnds - moveOnset;
            movSpks   = sum(spks >= moveOnset & spks <= moveEnds);
            frs(tr)   = movSpks / (movDur / 1000);
        end

        allFR(u, g) = mean(frs);
    end
end

% Normalize
peakFR     = max(allFR, [], 2);
validUnits = peakFR > 0;
allFRNorm  = allFR(validUnits, :) ./ peakFR(validUnits);
nValid     = sum(validUnits);

% Plotting
allFRNormClosed = [allFRNorm, allFRNorm(:, 1)];
popMean         = mean(allFRNorm, 1);
popMeanClosed   = [popMean, popMean(1)];

figure;
ax = polaraxes;
hold(ax, 'on');

% Individual neurons
for u = 1:nValid
    polarplot(ax, dirRadClosed, allFRNormClosed(u, :), ...
        'Color', [0.9 0.9 0.9], 'LineWidth', 0.3);
end

% Population mean
polarplot(ax, dirRadClosed, popMeanClosed, 'k-', 'LineWidth', 3);

title(ax, sprintf('Population Tuning Curves — Movement Period (n=%d neurons)', nValid));
hold(ax, 'off');
end