function [ciStruct, coverStruct, lengthStruct] = ...
    computeIntervalCoverageLength(...
        bsEstTableStruct, paramArray, targetParams, alphaCI)
% computeIntervalCoverageLength Computes bootstrap confidence intervals,
% coverage, and lengths.
%
% This function calculates confidence intervals for bootstrap estimates,
% checks if they cover the true parameter values, and computes the length
% of the intervals.
%
% Args:
%     bsEstTableStruct (struct): Struct containing tables of bootstrap 
%         debiased unit averaging estimates.
%         - Structure: paramName -> table of debiased estimates.
%         - Table columns correspond to different unit averaging approaches.
%     paramArray (cell array): Array of parameter definitions, where each 
%         struct includes:
%         - gradient (function handle): Function to compute gradients.
%         - saveName (string): Key used to save results for this parameter.
%     targetParams (struct): True parameter values. Structure: 
%          paramName -> true value.
%     alphaCI (double): Significance level for the confidence intervals. 
%          Confidence level is `1 - alphaCI`.
%
% Returns:
%     ciStruct (struct): Confidence intervals for target parameters based on
%         various averaging approaches. Structure: paramName -> 
%         approachName -> confidence interval.
%     coverStruct (struct): Boolean values indicating if intervals cover 
%         the true parameter value. Structure: paramName -> approachName 
%         -> Boolean.
%     lengthStruct (struct): Length of confidence intervals. 
%         Structure: paramName -> approachName -> length of interval.

    % Extract the number of parameters
    numParams = length(paramArray);

    % Initialize output structs
    ciStruct = struct();
    coverStruct = struct();
    lengthStruct = struct();

    % Loop through each parameter
    for paramID = 1:numParams
        % Extract parameter name
        paramName = paramArray{paramID}.saveName;

        % Extract true target value for the parameter
        targetValue = targetParams.(paramName);

        % Compute confidence intervals for all approaches
        confidenceIntervals = ...
            quantile(bsEstTableStruct.(paramName){:, :}, ...
                [alphaCI / 2, 1 - alphaCI / 2]);

        % Extract approach names
        approachNames = ...
            bsEstTableStruct.(paramName).Properties.VariableNames;

        % Loop through each averaging approach
        for approachID = 1:length(approachNames)
            % Get the name of the current approach
            approachName = string(approachNames{approachID});

            % Extract the current confidence interval
            currentCI = confidenceIntervals(:, approachID)';

            % Save the confidence interval to the output struct
            ciStruct.(paramName).(approachName) = currentCI;

            % Compute and save the length of the confidence interval
            lengthStruct.(paramName).(approachName) = ...
                currentCI(2) - currentCI(1);

            % Save whether the confidence interval covers the target value
            coverStruct.(paramName).(approachName) = ...
                (targetValue >= currentCI(1)) && ...
                (targetValue <= currentCI(2));
        end
    end
end