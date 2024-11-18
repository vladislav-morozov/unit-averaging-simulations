function outputWeights = ...
    weightsAICMMA(individualEstimators, y, x, weightScheme, numUnitsAvg)
    % weightsAICMMA Computes AIC and MMA weights for unit averaging.
    %
    % Calculates weights for averaging estimates from cross-sectional
    % units, using either Akaike Information Criterion (AIC) or Mallows
    % Model Averaging (MMA). The computed weights are designed to minimize
    % the prediction error for the target unit (first unit).
    %
    % Args:
    %     individualEstimators (matrix): k x N matrix of individual 
    %         estimates, where columns correspond to different 
    %         cross-sectional units.
    %     y (matrix): T x N matrix of outcome variables, with each column 
    %         representing the outcomes for a specific unit.
    %     x (array): T x N x 2 array of model covariates; the third 
    %         dimension indexes covariates.
    %     weightScheme (str): Scheme for weight calculation. Acceptable
    %         values are 'aic' (Akaike weights) or 'mma' (Mallows Model
    %         Averaging weights).
    %     numUnitsAvg (int, optional): Number of units to include in the 
    %         averaging calculation. Default is to use all units (N).
    %
    % Returns:
    %     outputWeights (vector): Weights vector of length k, computed 
    %         based on the specified weighting scheme.
    %
    % Example:
    %     weights = weightsAICMMA(individualEstimators, y, x, 'aic', 5);
    %     % Computes AIC weights using the first 5 units.
    
    % Extract dimensions from input
    [T, N, k] = size(x);

    % If numUnitsAvg is not provided, set it to N (use all units)
    if nargin < 5
        numUnitsAvg = N;
    end

    % Pre-allocate arrays for variance estimates and log-likelihoods
    sigmaHatSq = zeros(numUnitsAvg, 1);
    logLiks = zeros(numUnitsAvg, 1);

    % Covariates for the target unit (unit 1)
    covarsTarget = [x(:, 1, 1), x(:, 1, 2)];

    % Loop through the first `numUnitsAvg` units to calculate variances and
    % log-likelihoods
    for unitID = 1:numUnitsAvg
        % Covariates for the current unit
        covarsCurrent = [x(:, unitID, 1), x(:, unitID, 2)];

        % Compute residuals for current unit's estimate and observed values
        errorsVectorI = y(:, unitID) - ...
            covarsCurrent * individualEstimators(:, unitID);

        % Estimate variance for the current unit
        sigmaHatSq(unitID) = errorsVectorI' * errorsVectorI / (T - 2);

        % Residuals using current unit’s estimator applied to target unit's data
        errorsVector1 = y(:, 1) - ...
            covarsTarget * individualEstimators(:, unitID);

        % Compute log-likelihood based on the normal distribution
        logLiks(unitID) = -log(sigmaHatSq(unitID)) / 2 - ...
                          (1 / (2 * T)) * (errorsVector1' * errorsVector1)...
                          / sigmaHatSq(unitID);
    end

    % Calculate weights based on the selected weighting scheme
    if weightScheme == "aic"
        % AIC weights: Exponential transformation of the log-likelihoods
        outputWeights = exp(-2 * logLiks + 2 * k);
        % Normalize weights to sum to 1
        outputWeights = outputWeights / sum(outputWeights);

    else % MMA (Mallows Model Averaging) case
        % Construct matrix Z for MMA optimization (predictions for target unit)
        Z = nan(T, numUnitsAvg);
        for unitID = 1:numUnitsAvg
            Z(:, unitID) = covarsTarget * individualEstimators(:, unitID);
        end

        % Define covariance matrix D based on estimated variances
        D = diag(sigmaHatSq);

        % Set up options for quadratic programming
        options = optimoptions('quadprog', 'Display', 'none');

        % Solve the MMA weight optimization problem using quadprog
        outputWeights = quadprog(Z' * Z + D, -2 * y(:, 1)' * Z, [], [], ...
                                 ones(1, numUnitsAvg), 1, ...
                                 zeros(1, numUnitsAvg), [], [], options);
    end
end
