%% Part 1
%% Part A
nConditions = length(R);
nUnits  = length(R(1).unit);
startTime = -150;
endTime = 500;
binSize = 10;
bin_edges = (startTime:binSize:endTime);
nTimeBins = length(bin_edges) - 1; 
spikeCounts = zeros(nTimeBins, nConditions, nUnits);

for trial = 1:nTrials
    moveOnsetTime = R(trial).moveOnsetTime;
    for unit = 1:nUnits
        spikeTime = R(trial).unit(unit).spikeTimes;
        if isempty(spikeTime)
            continue
        end
        spikeTimeLocked = spikeTime - moveOnsetTime;
        spikeCounts(:, trial, unit) = histcounts(spikeTimeLocked, bin_edges);
    end
end

%% Part B
out = gaussFilt1(spikeCounts, 40);

%% Part C