function weightsRegArray = uaAverageWeight(paramArray, thetaPointArray, weightsArray, ...
thetaGrid, bandwidth) 
% uaAverageWeight Computes average weights for values in thetaGrid 

% Extract dimensions
numParams = length(paramArray);
numSamples = length(thetaPointArray);

% Approaches to plot
approachesToPlot = ["unrestr", "top",   "oracleSimilarity"];
 
% Loop through parameters
for paramID = 1:numParams
    % Extract parameter name
    paramName = paramArray{paramID}.saveName;
    
    % Combine together theta_1 samples drawn
    for sampleID = 1:numSamples
        if sampleID == 1
            thetasDrawn = thetaPointArray{sampleID}(1, :)';
        else
            thetasDrawn = [thetasDrawn; thetaPointArray{sampleID}(1, :)'];
        end
    end
    
    % Loop through the different approaches specified
    for approachID = 1:length(approachesToPlot)
        % Extract approach name
        approachName = approachesToPlot(approachID);
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
        weightsRegArray.(paramName).(approachName) = weightsNW;
    end
    
end

end
