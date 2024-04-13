function [errorStruct, estStruct, weightStruct, unitsUnrestrStruct] = ...
    uaSampleAveraging(y, covars,...
    thetaHat,thetaTrue,...
    estVarianceArray, ...
    paramArray, ...
    preparedMethodsArray, ...
    averagingMode)
% uaSampleAveraging Performs unit averaging and computes error. Applies
% averaging to all units that appear in the data, all parameters in
% paramArray. Applies all averaging methods specified by
% preparedMethodsArray
%
% Inputs: 1. y -- TxN matrix of outcome variables, columns index units
%         2. covars  -- TxNxk array of covariates that describe the model,
%            second dimension indexes units
%         3. thetaHat -- kxN matrix of estimated coefficients, columns
%            index units
%         4. thetaTrue -- kxN matrix of true coefficients, used to compute
%            true parameter values
%         5. estVarianceArray -- N-cell array where elements are covariance
%            matrices or kxkxN matrix, third index indices units
%         6. paramArray -- cell array of parameters. Each parameter must be
%               described by a struct that contains at least the .mu,
%               .gradient, and the .shortName fields
%         7. preparedMethodsArray -- cell array of averaging approaches.
%            Each approach is described by a struct that contains at least
%            the .shortName and the .unrestrictedArray fields
%         8. averagingMode -- 'all' or 'firstOnly'. Determines if all units
%            are to serve as targets or only the first one 
%
% Returns:
%         1. errorStruct -- struct with indexing structure .paramName ->
%            .approachName -> N-vector of errors. jth entry of vector is
%            error of estimating paramName with approach for unit j
%         2. estStruct
%         3. weightStruct -- struct with indexing structure .paramName ->
%            .approachName -> NxN matrix of weights. jth row of matrix
%            describes the weights when estimating paramName with
%            approachName for unit j.
%         4. unitsUnrestrStruct -- struct with indexing structure
%            .paramName -> .approachName -> NxN matrix of Booleans. jth
%            row of matrix describes units unrestricted when estimating
%            paramName with approachName for unit j

% Obtain dimensions
numParams = length(paramArray);
numApproaches = length(preparedMethodsArray);
if averagingMode == "all"
    numTargets = size(thetaHat, 2);
else
    numTargets = 1;
end
numUnits = size(thetaHat, 2);

% Preallocate suitable structs
for parID = 1:numParams
    % Extract parameter name
    paramName = paramArray{parID}.saveName;
    for approachID = 1:numApproaches
        % Extract averaging names
        approachShortName = preparedMethodsArray{approachID}.shortName;
        
        % Insert suitable nan arrays
        errorStruct.(paramName).(approachShortName) = nan(numTargets, 1);
        weightStruct.(paramName).(approachShortName) = ...
            nan(numTargets, numUnits);
        estStruct.(paramName).(approachShortName) = nan(numTargets, 1);
    end
end


% Averaging: loop through parameters, targets, approaches
for parID = 1:numParams
    
    % Extract current target parameter
    paramStruct = paramArray{parID};
    % Extract parameter name
    paramName = paramArray{parID}.saveName;
    % Construct all individual estimates
    unitParamEsts = paramStruct.mu(thetaHat, covars, y)';
    
    % Loop through target units
    for targetID = 1:numTargets
        
        % Compute target value of parameter for the target
        targetParamValue = paramStruct.mu(...
            thetaTrue(:, targetID), ...
            covars(:, targetID,:), y(:, targetID));
        
        % Estimate the gradient for the target unit
        gradientEst = paramStruct.gradient(...
            thetaTrue(:, targetID), ...
            covars(:, targetID,:), y(:, targetID));
        
        % Loop through averaging approaches
        for approachID = 1:numApproaches
            % Extract name
            approach = preparedMethodsArray{approachID};
            approachShortName = approach.shortName;
            
            % Create default vector of unrestricted units
            unrestrictedBooltargetId = true(numTargets, 1);
            
            % if approach is optimal, extract its restricted component and
            % modify the unrestrictedBooltargetId vector
            if isfield(approach, "unrestrictedArray")
                % If the approach is 'top', then further differentiate
                if approach.shortName(1:3) == "top" 
                % Top approaches are based on fixed-N weights
                    targetWeightsSetup ...
                        = weightStruct.(paramName).unrestr(targetID, :)';
                else
                % Other approaches get a dummy weight
                    targetWeightsSetup = ones(numTargets, 1)/numTargets;
                end
                
                % Extract the unrestricted units
                unrestrictedBooltargetId = ...
                    approach.unrestrictedArray{targetID}(targetWeightsSetup);
                % Save the unrestricted units
                unitsUnrestrStruct.(paramName).(approachShortName)(...
                    targetID, :) = ...
                    unrestrictedBooltargetId';
            end
           
            % Compute weights of the current approach
            weightsCurrent = approach.weightFunction(...
                y, covars, thetaHat,...
                estVarianceArray, gradientEst, ...
                targetID, unrestrictedBooltargetId);
            
            % Compute averaging estimate and its error
            avgEst = weightsCurrent'*unitParamEsts;  
            errorCurrent = avgEst-targetParamValue;
            
            % Insert the weights into the weights array
            weightStruct.(paramName).(approachShortName)(targetID, :) = ...
                weightsCurrent';
            % Insert the error into the error array
            errorStruct.(paramName).(approachShortName)(targetID) = ...
                errorCurrent;
            % Insert the estimate into the estimate array
            estStruct.(paramName).(approachShortName)(targetID) = avgEst;
        end
    end
end
end