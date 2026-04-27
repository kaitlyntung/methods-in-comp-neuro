%% Part A
num_trials = length(Kinematics.Trials);

for trial = 1:num_trials
    Kinematics.Trials(trial) = zscore(Kinematics.Trials(trial));
end