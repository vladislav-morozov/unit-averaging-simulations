function [errorStruct, estStruct, weightStruct, unitsUnrestrStruct] = ...
    sampleAveraging(y, covars, ...
        thetaHat, thetaTrue, ...
        estVarianceArray, ...
        paramArray, ...
        preparedMethodsArray, ...
        averagingMode)
% sampleAveraging Performs unit averaging, computes errors, and returns
% weights and unrestricted units.
%
% This function calculates weighted averages of parameter estimates for
% target units using various averaging methods. It computes estimation
% errors, stores the weights used, and tracks which units were unrestricted
% during averaging.
%
% Args:
%     y (TxN matrix): Outcome variables. Rows index time periods, columns 
%         index units.
%     covars (TxNxk array): Covariates describing the model. 
%         - Rows index time periods.
%         - Columns index units.
%         - Third dimension indexes covariate variables.
%     thetaHat (kxN matrix): Estimated coefficients. 
%         - Rows index coefficient dimensions.
%         - Columns index units.
%     thetaTrue (kxN matrix): True coefficients used to compute actual 
%         parameter values.
%     estVarianceArray (cell array or kxkxN matrix): 
%         - Covariance matrices for each unit. Can be a cell array 
%         (N elements) or a 3D matrix.
%     paramArray (cell array): Array of parameter structs, where each 
%         struct includes:
%         - mu (function handle): Function to compute parameter estimates.
%         - gradient (function handle): Function to compute parameter 
%         gradients.
%         - shortName (string): Short identifier for the parameter.
%         - saveName (string): Key used to save results for this parameter.
%     preparedMethodsArray (cell array): Array of averaging method structs,
%         where each struct includes:
%         - shortName (string): Short identifier for the method.
%         - unrestrictedArray (cell array, optional): Determines 
%         unrestricted units for the method.
%         - weightFunction (function handle): Function to compute averaging
%         weights.
%     averagingMode (string): Determines target units for averaging:
%         - 'all': All units serve as targets.
%         - 'firstOnly': Only the first unit is a target.
%
% Returns:
%     errorStruct (struct): Nested structure indexed as `paramName -> 
%       approachName -> J-vector of errors`. Here J = 1 is averagingMode is
%       equal to 'firstOnly'; J = N if averagingMode = 'all'.
%         - Each vector's j-th entry corresponds to the estimation error 
%        for target unit j.
%     estStruct (struct): Nested structure indexed as `paramName -> 
%         approachName -> k-vector of estimates`.
%         - Each vector's j-th entry corresponds to the parameter estimate 
%         for target unit j.
%     weightStruct (struct): Nested structure indexed as `paramName -> 
%         approachName -> k x N matrix of weights`.
%         - Each row describes weights applied to each unit for a specific 
%         target.
%     unitsUnrestrStruct (struct): Nested structure indexed as `paramName 
%         -> approachName -> k x N matrix of Booleans`.
%         - Each row describes unrestricted units for a specific target.

    % Obtain dimensions
    numParams = length(paramArray);
    numApproaches = length(preparedMethodsArray);
    numUnits = size(thetaHat, 2);

    if averagingMode == "all"
        numTargets = numUnits;
    else
        numTargets = 1;
    end

    % Preallocate structures for results
    for parID = 1:numParams
        % Extract parameter name
        paramName = paramArray{parID}.saveName;
        for approachID = 1:numApproaches
            % Extract averaging approach name
            approachShortName = preparedMethodsArray{approachID}.shortName;

            % Initialize result arrays with NaNs
            errorStruct.(paramName).(approachShortName) = ...
                nan(numTargets, 1);
            weightStruct.(paramName).(approachShortName) = ...
                nan(numTargets, numUnits);
            estStruct.(paramName).(approachShortName) = nan(numTargets, 1);
        end
    end

    % Perform averaging
    for parID = 1:numParams
        % Extract parameter details
        paramStruct = paramArray{parID};
        paramName = paramStruct.saveName;

        % Compute individual unit estimates
        unitParamEsts = paramStruct.mu(thetaHat, covars, y)';

        % Extract the parameter function
        paramFun = @(theta) paramStruct.mu(theta, covars, y);

        % Loop through target units
        for targetID = 1:numTargets
            % Compute target parameter value
            targetParamValue = paramStruct.mu(...
                thetaTrue(:, targetID), ...
                covars(:, targetID, :), ...
                y(:, targetID));

            % Compute parameter gradient for the target unit
            gradientEst = paramStruct.gradient(...
                thetaTrue(:, targetID), ...
                covars(:, targetID, :), ...
                y(:, targetID));

            % Loop through averaging approaches
            for approachID = 1:numApproaches
                % Extract averaging method details
                approach = preparedMethodsArray{approachID};
                approachShortName = approach.shortName;

                % Initialize unrestricted units as true by default
                unrestrictedBool = true(numUnits, 1);

                % Handle optimal methods with unrestricted units
                if isfield(approach, "unrestrictedArray")
                    if approachShortName(1:3) == "top"
                        % Top methods: Use fixed weights
                        targetWeightsSetup = ...
                            weightStruct.(paramName).unrestr(targetID, :)';
                    else
                        % Default weights for non-top methods
                        targetWeightsSetup = ones(numUnits, 1) / numUnits;
                    end

                    % Determine unrestricted units
                    unrestrictedBool = ...
                        approach.unrestrictedArray{targetID}(...
                        targetWeightsSetup, paramFun);

                    % Store unrestricted units
                    unitsUnrestrStruct.(paramName).(approachShortName)(...
                        targetID, :) = ...
                        unrestrictedBool';
                end

                % Compute weights for the current approach
                weightsCurrent = approach.weightFunction(...
                    y, covars, thetaHat, ...
                    estVarianceArray, gradientEst, ...
                    targetID, unrestrictedBool);

                % Compute averaged estimate and error
                avgEst = weightsCurrent' * unitParamEsts;
                errorCurrent = avgEst - targetParamValue;

                % Store results
                weightStruct.(paramName).(approachShortName)(targetID, :) = ...
                    weightsCurrent';
                errorStruct.(paramName).(approachShortName)(targetID) = ...
                    errorCurrent;
                estStruct.(paramName).(approachShortName)(targetID) = avgEst;
            end
        end
    end
end
