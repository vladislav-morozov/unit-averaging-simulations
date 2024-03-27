function [weightsRegStruct, averageFirstWeightStruct] ...
    = uaAverageWeight(paramArray, thetaPointArray, weightsArray, ...
thetaGrid, bandwidth) 
% uaAverageWeight Computes average weights for values in thetaGrid. Also
% computes average own weight

% Extract dimensions
numParams = length(paramArray);
numSamples = length(thetaPointArray);

% Approaches to plot
approachesAverageWeights = ["unrestr", "top",   "oracleSimilarity"];
% Approaches for own weight
approachesOwnWeight = ["unrestr", "top", "oracleSimilarity", "stein"];
 
% Loop through parameters
for paramID = 1:numParams
    % Extract parameter name
    paramName = paramArray{paramID}.saveName;
    
    %%%%%%%%%%%%%%%%%%%%%%%
    %%% Average weights %%%
    %%%%%%%%%%%%%%%%%%%%%%%
    
    % Combine together theta_1 samples drawn
    for sampleID = 1:numSamples
        if sampleID == 1
            thetasDrawn = thetaPointArray{sampleID}(1, :)';
        else
            thetasDrawn = [thetasDrawn; thetaPointArray{sampleID}(1, :)'];
        end
    end
    
    % Loop through the different approaches specified
    for approachID = 1:length(approachesAverageWeights)
        % Extract approach name
        approachName = approachesAverageWeights(approachID);
        % Combine together the weights obtain (same order as thetasDrawn)
        for sampleID = 1:numSamples
            if sampleID == 1
                weightsDrawn = weightsArray{sampleID}.(paramName).(approachName)';
            else
                weightsDrawn = [weightsDrawn; weightsArray{sampleID}.(paramName).(approachName)'];
            end
        end
        
        % Run Nadaraya-Watson regressions of weights on points
        for pointID = 1:length(thetaGrid)
            kernelWeights = exp(-(thetasDrawn - thetaGrid(pointID)).^2 / (2 * bandwidth^2));
            weightsNW(pointID) = sum(weightsDrawn .* kernelWeights) / sum(kernelWeights);
        end
        % Save computed weights into a suitable struct
        weightsRegStruct.(paramName).(approachName) = weightsNW;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Average first weight %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for approachID = 1:length(approachesOwnWeight)
        % Extract approach name
        approachName = approachesOwnWeight(approachID);
        % Extract weight of the first unit in each iteration
        for sampleID = 1:numSamples
            if sampleID == 1
                firstWeight = weightsArray{sampleID}.(paramName).(approachName)(1);
            else
                firstWeight = [firstWeight; weightsArray{sampleID}.(paramName).(approachName)(1)];
            end
        end
        % Save average first unit weight for current approach
        averageFirstWeightStruct.(paramName).(approachName) = mean(firstWeight);
    end
end

end
