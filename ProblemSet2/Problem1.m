%% Part A
neurons = length(Data.spikes);
mean_spikes = cell(1, neurons);

for neuron = 1:neurons
    trial_mask = Data.stimuli{1, neuron}(:, 2) == 700;
    spikes = Data.spikes{1, neuron}(trial_mask);
    stimuli = Data.stimuli{1, neuron}(trial_mask, 3);
    unique_conds = unique(stimuli);
    for cond = 1:length(unique_conds)
        cond_mask = stimuli == unique_conds(cond);
        spikes_per_cond = spikes(cond_mask);
        counts = cellfun(@length, spikes_per_cond);
        mean_spikes{neuron}(cond) = mean(counts);
    end
end

%% Part B
