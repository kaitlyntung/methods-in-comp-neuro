function [firingRateMatrix, conditionLabels] = buildFiringRateMatrix(R)
% Build a nTrials x nUnits matrix of mean firing rates during the
% full movement period for each trial, z-scored across trials. Inputs: maze
% data struc

nTrials = length(R);
nUnits = length(R(1).unit);

% Compute reach direction group for each trial
targetXY = reshape([R.targetXY], 2, [])';
angles = atan2d(targetXY(:,2), targetXY(:,1));
angleGroups = round(angles / 45) * 45;
angleGroups(angleGroups == 180) = -180;
angleGroups(angleGroups == 0 & angles < 5)  = -45;
angleGroups(angleGroups == 0 & angles >= 5) = 45;

% Pre-allocate
firingRateMatrix = zeros(nTrials, nUnits);

for tr = 1:nTrials
    moveOnset = R(tr).moveOnsetTime;
    moveEnds = R(tr).moveEndsTime;
    movDur = moveEnds - moveOnset;

    for u = 1:nUnits
        spks = R(tr).unit(u).spikeTimes;
        alignedSpks = spks - moveOnset;

        % Count spikes over full movement period
        movSpks = sum(alignedSpks >= 0 & alignedSpks <= movDur);
        firingRateMatrix(tr, u) = movSpks / (movDur / 1000);
    end
end

% Z-score each neuron across trials
firingRateMatrix = zscore(firingRateMatrix);

conditionLabels = angleGroups;
end