function [averageWeightsGridStruct, ...
    averageFirstWeightStruct,...
    averageMaxDiff, ...
    averageMassDiff, ...
    averageUnitsUnrestrStruct] ...
    = processOneGridValueWeightsUnrestUnits(paramArray, thetaPointArray,...
    weightsArray, unitsUnrestrArray, ...
    thetaGrid, bandwidth) 
% processOneGridValueWeightsUnrestUnits Processes Monte Carlo results for
% weights and unrestricted units at a specific grid point.
%
% This function estimates expected weights, probabilities of being
% unrestricted, and differences in weight and mass between restricted and
% unrestricted units across Monte Carlo samples for various grid points in
% the parameter space.
%
% Args:
%     paramArray (cell array): Array of structs, each describing a 
%         parameter. Must contain .saveName (string): Key used to save 
%         results for this parameter.
%     thetaPointArray (cell array): True values of thetas for each Monte 
%         Carlo sample. Each cell contains a 2xN array of theta values for
%         a specific sample.
%     weightsArray (cell array): Estimated weights for each sample. Each 
%         cell is a struct with the structure: paramName -> approachName 
%         -> weight vector.
%     unitsUnrestrArray (cell array): Indicators of unrestricted units for
%          each sample. Same structure as `weightsArray`.
%     thetaGrid (vector): Grid points at which results are reported.
%     bandwidth (scalar): Bandwidth used in Nadaraya-Watson regression.
%
% Returns:
%     averageWeightsGridStruct (struct): Expected weights for units with a 
%         given parameter value. 
%         - Structure: paramName -> approachName -> weight vector 
%         (length of `thetaGrid`).
%     averageFirstWeightStruct (struct): Expected weight for the target 
%         unit itself.
%         - Structure: paramName -> approachName -> scalar (average weight).
%     averageMaxDiff (struct): Expected maximum absolute difference between 
%         unrestricted and restricted weights.
%         - Structure: paramName -> approachName -> scalar (average max
%         difference).
%     averageMassDiff (struct): Average difference in total mass between 
%         unrestricted and restricted units.
%         - Structure: paramName -> approachName -> scalar (average mass 
%         difference).
%     averageUnitsUnrestrStruct (struct): Probability a unit is  
%         unrestricted for a given parameter value.
%         - Structure: paramName -> approachName -> vector (length of 
%         `thetaGrid`).


% Extract dimensions
numParams = length(paramArray);
numSamples = length(thetaPointArray);
 
% Loop through parameters
for paramID = 1:numParams
    % Extract parameter name
    paramName = paramArray{paramID}.saveName;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%          Expected weights         %%%
    %%% Probability of being unrestricted %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Extract names of approaches present
    approachesAverageWeights = fields(unitsUnrestrArray{1}.(paramName));
    approachesOwnWeight = approachesAverageWeights;
    
    % Combine together theta_1 samples drawn into a single array
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
        approachName = approachesAverageWeights{approachID};
        % Combine together the weights obtain (same order as thetasDrawn)
        for sampleID = 1:numSamples
            if sampleID == 1
                weightsDrawn = ...
                    weightsArray{sampleID}.(paramName).(approachName)';
                unitsUnrestr = ...
                    unitsUnrestrArray{sampleID}.(paramName).(approachName)';
            else
                weightsDrawn = ...
                    [weightsDrawn; ...
                      weightsArray{sampleID}.(paramName).(approachName)'];
                unitsUnrestr = ...
                    [unitsUnrestr; ...
                    unitsUnrestrArray{sampleID}.(paramName).(approachName)'];
            end
        end
        
        % Run Nadaraya-Watson regressions of weights on grid points, this
        % estimates the expected weights as a function of underlying value
        % of thetas. Same argument estimates probability of being
        % unrestricted as a function of parameter value
        for pointID = 1:length(thetaGrid)
            kernelWeights = exp(-(thetasDrawn - thetaGrid(pointID)).^2 / (2 * bandwidth^2));
            weightsNW(pointID) = sum(weightsDrawn .* kernelWeights) / sum(kernelWeights);
            unitsNW(pointID) = sum(unitsUnrestr .* kernelWeights) / sum(kernelWeights);
        end
        % Save computed weights into a suitable struct
        averageWeightsGridStruct.(paramName).(approachName) = weightsNW;
        averageUnitsUnrestrStruct.(paramName).(approachName) = unitsNW;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Average first weight %%%
    %%% Max diff in weights  %%%
    %%% Difference in mass   %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Loop through approaches
    for approachID = 1:length(approachesOwnWeight)
        % Extract approach name
        approachName = approachesOwnWeight{approachID};
        % Extract weight of the first unit in each iteration
        for sampleID = 1:numSamples
            if sampleID == 1
                % Extract current weight vector and restricted units
                currentWeights = ...
                    weightsArray{sampleID}.(paramName).(approachName);
                currentRestr = ...
                    ~unitsUnrestrArray{sampleID}.(paramName).(approachName);
                
                % Extract weights for the fixed-N approach
                currentWeightsFixedN = ...
                    weightsArray{sampleID}.(paramName).unrestr;
                
                % Extract weight of unit 1
                firstWeight = currentWeights(1);
                
                % Compute max absolute difference in restricted units
                restrWeightDiff = currentWeights(currentRestr)-...
                    currentWeightsFixedN(currentRestr);
                
                % Extract maximum difference
                maxDiff = max(abs(restrWeightDiff));
                
                % Exract difference in mass
                massDiff = sum(restrWeightDiff);
            else
                % Extract current weight vector and restricted units
                currentWeights = ...
                    weightsArray{sampleID}.(paramName).(approachName);
                currentRestr = ...
                    ~unitsUnrestrArray{sampleID}.(paramName).(approachName);
                
                % Extract weights for the fixed-N approach
                currentWeightsFixedN = ...
                    weightsArray{sampleID}.(paramName).unrestr;
               
                % Compute max absolute difference in restricted units
                restrWeightDiff = currentWeights(currentRestr)-...
                    currentWeightsFixedN(currentRestr); 
                
                % Add to arrays: first weights, max and mass differences
                firstWeight = [firstWeight; currentWeights(1);]; 
                maxDiff = [maxDiff; max(abs(restrWeightDiff))];
                massDiff = [massDiff; sum(restrWeightDiff);];
            end
        end
        % Save average first unit weight for current approach
        averageFirstWeightStruct.(paramName).(approachName) = ...
            mean(firstWeight);
        % Save max differences 
        averageMaxDiff.(paramName).(approachName) = mean(maxDiff);
        % Save differences in masses
        averageMassDiff.(paramName).(approachName) = mean(massDiff);
    end
end
end
