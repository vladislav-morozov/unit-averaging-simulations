function estDebiasedStructBs = ...
    debiasedBootstrapEstimator(...
        estStructBs, weightStructBs, ...
        paramArray, ...
        thetaHatBs, covarsBs, yBs ...
       )
% debiasedBootstrapEstimator Computes plug-in debiased unit averaging
% estimators for supplied bootstrap sample.
%
% This function adjusts unit-averaged parameter estimates to reduce bias,
% based on a provided bootstrap sample. The debiasing uses gradient-based
% bias correction with respect to the target unit (assumed to be the first
% unit).
%
% Args:
%     estStructBs (struct): Estimates of parameters for unit averaging 
%       approaches. Structure: paramName -> averagingApproach -> 
%       estimate value.
%     weightStructBs (struct): Weights for parameters under different unit
%       averaging approaches. Structure: paramName -> averagingApproach 
%       -> weight vector.
%     paramArray (cell array): Array of parameter structs, where each 
%        struct includes:
%         - gradient (function handle): Function to compute gradients.
%         - saveName (string): Key used to save results for this parameter.
%     thetaHatBs (matrix): k x N  matrix of estimated coefficients for the 
%        bootstrap sample. Columns correspond to individual units.
%     covarsBs (array): T x N x k array of covariates for 
%        the bootstrap sample. Second dimension indexes units.
%     yBs (matrix): T x N matrix of outcomes for the bootstrap sample.
%         Columns correspond to individual units.
%
% Returns:
%     estDebiasedStructBs (struct): Debiased unit averaging estimators for 
%         the current bootstrap sample. Structure: paramName -> 
%         averagingApproach -> debiased estimate.

    % Extract the number of parameters and the sample size
    numParams = length(paramArray);
    T = size(yBs, 1);

    % Initialize output struct
    estDebiasedStructBs = struct();

    % Loop through parameters
    for paramID = 1:numParams
        % Extract parameter name
        paramName = paramArray{paramID}.saveName;

        % Get the list of unit averaging approaches for this parameter
        approaches = fieldnames(estStructBs.(paramName));

        % Loop through averaging approaches
        for approachID = 1:length(approaches)
            % Extract approach name
            approachName = string(approaches{approachID});

            % Retrieve the weights for the current approach
            weights = weightStructBs.(paramName).(approachName);

            % Compute the gradient for the target unit (assumed to be the
            % first unit)
            gradientTarget = paramArray{paramID}.gradient(...
                thetaHatBs(:, 1), covarsBs(:, 1, :), yBs(:, 1));

            % Compute the bias correction term
            biasCorrection = ...
                weights * (gradientTarget' * ...
                (thetaHatBs - thetaHatBs(:, 1)))' / sqrt(T);

            % Apply the bias correction to the estimate
            estDebiasedStructBs.(paramName).(approachName) = ...
                estStructBs.(paramName).(approachName) - biasCorrection;
        end
    end
end
