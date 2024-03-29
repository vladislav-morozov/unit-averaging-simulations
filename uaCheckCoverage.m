function [ciStruct, coverStruct, lengthStruct] = uaCheckCoverage(...
    bsEstTableStruct, paramArray, targetParams, alphaCI)
% uaCheckCoverage Computes CIs and checks if the cover the target value

% Extract number of parameters
numParams = length(paramArray);

% Loop through parameters
for paramID = 1:numParams
    % Extract parameter name
    paramName = paramArray{paramID}.saveName;
    
    % Target value
    targetValue = targetParams.(paramName);
    
    % Compute confidence intervals
    confidenceIntervals = ...
        quantile(bsEstTableStruct.(paramName){:, :}, ...
        [alphaCI/2, 1-alphaCI/2]);
    
    % Approach names
    approachNames = bsEstTableStruct.(paramName).Properties.VariableNames;
     
    % Loop through approaches
    for approachID = 1:length(approachNames)
       % Approach name
       approachName = string(approachNames{approachID});
       
       % Current CI
       currentCI = confidenceIntervals(:, approachID)';
       
       % Save CI
       ciStruct.(paramName).(approachName) = currentCI;
       
       % Compute interval length
       lengthStruct.(paramName).(approachName) = currentCI(2)-currentCI(1);
       
       % Compute whether the interval contains the target value
       coverStruct.(paramName).(approachName) = ...
        (targetValue>=currentCI(1)) && ...
             (targetValue<=currentCI(2)); 
    end
    
    
end
end