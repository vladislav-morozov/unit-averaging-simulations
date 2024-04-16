function [weightsRegStruct, averageFirstWeightStruct,...
    averageMaxDiff, averageMassDiff, ...
    unitsUnrestrStruct] ...
    = uaAverageDiffWeightsUnrestrictedUnits(paramArray, thetaPointArray,...
    weightsArray, unitsUnrestrArray, ...
    thetaGrid, bandwidth) 
% uaAverageWeight Computes average weights for values in thetaGrid. Also
% computes average own weight

% Extract dimensions
numParams = length(paramArray);
numSamples = length(thetaPointArray);

% Approaches to plot
approachesAverageWeights = ...
    ["unrestr", ...
     "focus_oracle_10", "focus_oracle_25", ...
     "focus_oracle_10_pct","focus_oracle_25_pct", "focus_oracle_50_pct"...
     "focus_cluster_2", "focus_cluster_4", "focus_cluster_8", ...
     "stein", ...
     "top_10", "top_25", "top_10_pct", "top_25_pct", "top_50_pct",...
     "random_10", "random_20", ...
     "oracleClasses", "antiOracleClasses", ...
      "oracle_similar_10", "oracle_similar_25" , ...
     "oracle_similar_10_pct", "oracle_similar_25_pct", ...
     "oracle_similar_50_pct", ... 
     "cluster_coef_2", "cluster_coef_4", "cluster_coef_8", ...
   ];

% Approaches for own weight
approachesOwnWeight = approachesAverageWeights;
 
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
                unitsUnrestr = unitsUnrestrArray{sampleID}.(paramName).(approachName)';
            else
                weightsDrawn = [weightsDrawn; weightsArray{sampleID}.(paramName).(approachName)'];
                unitsUnrestr = [unitsUnrestr; unitsUnrestrArray{sampleID}.(paramName).(approachName)'];
            end
        end
        
        % Run Nadaraya-Watson regressions of weights on points
        for pointID = 1:length(thetaGrid)
            kernelWeights = exp(-(thetasDrawn - thetaGrid(pointID)).^2 / (2 * bandwidth^2));
            weightsNW(pointID) = sum(weightsDrawn .* kernelWeights) / sum(kernelWeights);
            unitsNW(pointID) = sum(unitsUnrestr .* kernelWeights) / sum(kernelWeights);
        end
        % Save computed weights into a suitable struct
        weightsRegStruct.(paramName).(approachName) = weightsNW;
        unitsUnrestrStruct.(paramName).(approachName) = unitsNW;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Average first weight %%%
    %%% Max diff in weights  %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for approachID = 1:length(approachesOwnWeight)
        % Extract approach name
        approachName = approachesOwnWeight(approachID);
        % Extract weight of the first unit in each iteration
        for sampleID = 1:numSamples
            if sampleID == 1
                % Extract current weight vector and restricted units
                currentWeights = weightsArray{sampleID}.(paramName).(approachName);
                currentRestr = ~unitsUnrestrArray{sampleID}.(paramName).(approachName);
                
                % Extract weights for the fixed-N approach
                currentWeightsFixedN = weightsArray{sampleID}.(paramName).unrestr;
                
                % Extract weight of unit 1
                firstWeight = currentWeights(1);
                
                % Compute max absolute difference in restricted units
                restrWeightDiff = currentWeights(currentRestr)-...
                    currentWeightsFixedN(currentRestr);
                
                maxDiff = max(abs(restrWeightDiff));
                massDiff = sum(restrWeightDiff);
            else
                % Extract current weight vector and restricted units
                currentWeights = weightsArray{sampleID}.(paramName).(approachName);
                currentRestr = ~unitsUnrestrArray{sampleID}.(paramName).(approachName);
                
                % Extract weights for the fixed-N approach
                currentWeightsFixedN = weightsArray{sampleID}.(paramName).unrestr;
                
                 
                
                % Compute max absolute difference in restricted units
                restrWeightDiff = currentWeights(currentRestr)-...
                    currentWeightsFixedN(currentRestr); 
                
                % Add to arrays
                firstWeight = [firstWeight; currentWeights(1);]; %#ok<*AGROW>
                maxDiff = [maxDiff; max(abs(restrWeightDiff))];
                massDiff = [massDiff; sum(restrWeightDiff);];
            end
        end
        % Save average first unit weight for current approach
        averageFirstWeightStruct.(paramName).(approachName) = mean(firstWeight);
        averageMaxDiff.(paramName).(approachName) = mean(maxDiff);
        averageMassDiff.(paramName).(approachName) = mean(massDiff);
    
        
    
    end
    
    
end

end
