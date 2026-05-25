function decoderModel = trainDecoder(R, nFolds)
% Bin neural data and computes hand velocity
% for all trials, then trains ridge and random forest decoders.
% Returns both trained models in a single struc

nTrials   = length(R);
nUnits    = length(R(1).unit);
binSizeMs = 100;

neuralData   = cell(nTrials, 1);
velocityData = cell(nTrials, 1);

for tr = 1:nTrials
    moveOnset = R(tr).moveOnsetTime;
    moveEnds  = R(tr).moveEndsTime;
    binEdges  = moveOnset : binSizeMs : moveEnds;
    nBins     = length(binEdges) - 1;

    if nBins < 1 
        continue 
    end

    % Bin spike counts
    trialNeural = zeros(nBins, nUnits);
    for u = 1:nUnits
        spks = R(tr).unit(u).spikeTimes;
        trialNeural(:, u) = histcounts(spks, binEdges) / (binSizeMs / 1000);
    end

    trialNeural = sqrt(trialNeural);

    % Compute hand velocity
    handTimes  = R(tr).HAND.times;
    dt         = mean(diff(handTimes));
    velX       = gradient(R(tr).HAND.X, dt);
    velY       = gradient(R(tr).HAND.Y, dt);
    binCenters = binEdges(1:end-1) + binSizeMs / 2;

    trialVel      = zeros(nBins, 2);
    trialVel(:,1) = interp1(handTimes, velX, binCenters);
    trialVel(:,2) = interp1(handTimes, velY, binCenters);

    neuralData{tr}   = trialNeural;
    velocityData{tr} = trialVel;
end

% Train both decoders
validTrials  = ~cellfun(@isempty, neuralData);
neuralData   = neuralData(validTrials);
velocityData = velocityData(validTrials);

ridgeModel = trainRidgeDecoder(neuralData, velocityData, nFolds);
rfModel    = trainRFDecoder(neuralData, velocityData, nFolds);

decoderModel.ridge     = ridgeModel;
decoderModel.rf        = rfModel;
decoderModel.binSizeMs = binSizeMs;
end