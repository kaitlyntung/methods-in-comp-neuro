%% Part A
start_time = 100;
end_time = 400;


num_trials = length(R);
num_units = length(R(1).unit);
spikeCounts = zeros(num_units, num_trials);

for trial = 1:num_trials
    targetAppearsTime = R(trial).targetAppearsTime;
    start_window = targetAppearsTime + start_time;
    end_window = targetAppearsTime + end_time;
    for unit = 1:num_units
        spike = R(trial).unit(unit).spikeTimes;
        spikeCounts(unit, trial) = sum(spike > start_window & spike <= end_window);
    end
end
%% Part B
window_size = (abs(end_time) + abs(start_time)) * 0.001;
conditions = [R.conditionID];
cond_one = conditions(1);
cond_two = conditions(2);
cond_one_idx = R.conditionID == cond_one;
cond_two_idx = R.conditionID == cond_two;

avgRateOne = mean(spikeCounts(:, cond_one_idx), 2) / window_size;
avgRateTwo = mean(spikeCounts(:, cond_two_idx), 2) / window_size;

%% Part C
diff       = abs(meanRate1 - meanRate2);
tunedMask  = diff >= 2.0;
tunedIdx   = find(tunedMask);

fprintf('\nPart 1C — Tuned units (|ΔFR| ≥ 2 sp/s): %d / %d\n', ...
        sum(tunedMask), nUnits);
fprintf('\n%5s  %14s  %14s  %8s\n', 'Unit', ...
        sprintf('Cond%d (sp/s)', cond1), ...
        sprintf('Cond%d (sp/s)', cond2), '|Diff|');
fprintf('%s\n', repmat('-', 1, 48));
for k = 1:length(tunedIdx)
    u = tunedIdx(k);
    fprintf('%5d  %14.2f  %14.2f  %8.2f\n', u, meanRate1(u), meanRate2(u), diff(u));
end