function estDebiasedStructBs = ...
    uaDebiasedBootstrapEvaluator(...
        estStructBs, weightStructBs, ...
        paramArray, ...
        thetaHatBs, covarsBs, yBs ...
       )
% Data sizes
numParams = length(paramArray);
T = size(yBs, 1);

for paramID = 1:numParams
    % Extract parameter name
    paramName = paramArray{paramID}.saveName ;
    
    % Extract approaches
    approaches = fieldnames(estStructBs.(paramName));
    
    % Loop through approaches
    for approachID = 1:length(approaches)
        % Extract approach name
        approachName = string(approaches{approachID});
        
        % Current weights
        weights = weightStructBs.(paramName).(approachName);
        
        % Current gradient
        d1 = ...
            paramArray{paramID}.gradient(...
            thetaHatBs(:, 1), covarsBs(:, 1, :), yBs(:, 1)...
            );
        
        % Bias correction term
        biasCorrection = ...
            weights*(d1'*(thetaHatBs-thetaHatBs(:, 1)))'/sqrt(T);
        
        % Corrected estimate
        estDebiasedStructBs.(paramName).(approachName) = ...
            estStructBs.(paramName).(approachName) - biasCorrection;
    end
end

end