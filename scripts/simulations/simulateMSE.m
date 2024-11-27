% ===========================================================
% File: simulateMSE.m
% Description: This script simulates the Mean Squared Error (MSE),
%              bias, variance, weight, and unrestricted unit properties
%              of unit averaging approaches in a heterogeneous 
%              panel ARX(1) model. 
%
% ===========================================================
%
% Project Name: Unit Averaging for Heterogeneous Panels
% Authors: Christian Brownlees, Vladislav Morozov
%
% Model: Heterogeneous panel ARX(1) model:
%   y_{it} = theta_{i1} y_{it-1} + theta_{i2} x_{it} + u_{it}
%
% Simulation Workflow:
% 1. For each combination of N and T (number of units, time periods):
%    a. Iterate over candidate values of theta_{i1} (grid points).
%      I. Iterate over samples for each candidate value:
%           i. Generate data and estimate parameters.
%           ii. Apply unit averaging schemes to estimate target parameters.
%           iii. Record averaging errors and predictions.
%           iv. Record averaging weights and unrestricted units.
%      II. Aggregate results for the grid point (e.g., compute MSE).
%    b. Aggregate results over grid into a common format (see below).
% 2. Save and export results.
% 
% Output format:
% - Results are saved in cell arrays indexed by combinations of (N, T).
% - Results come in two formats:
%     - Tables with columns indexing unit averaging approaches. Applies to 
%       MSE, bias, variance, differences in masses between unrestricted and
%       restricted units, average own weight.
%     - Structs indexed by averaging approaches. Applies to average weights
%       and probability of being unrestricted.
% ===========================================================

%% Simulation Parameters and Initialization

% Extract the number of dimensions for each parameter
numT = length(valuesT);                % Number of time periods
numN = length(valuesN);                % Number of units
numParams = length(paramArray);        % Number of parameters
numTargetPoints = length(theta1Range); % Grid points for theta_{i1}

% Preallocate result storage for each (N, T) combination
mseTablesNT = cell(numN, numT);        % MSE results
biasTablesNT = cell(numN, numT);       % Bias results
varTablesNT = cell(numN, numT);        % Variance results
weightsTablesNT = cell(numN, numT);    % All weights
firstWeightNT = cell(numN, numT);      % Weights for the first unit
maxDiffNT = cell(numN, numT);          % Maximum difference in weights
massDiffNT = cell(numN, numT);         % Mass difference (restr vs unrestr)
unitsUnrestrNT = cell(numN, numT);     % Probability of being unrestricted

%% Main simulation loop
 
% Loop over the number of time periods (T) and number of units (N)
for tID = 1:numT
    for nID=1:numN
        % Current sample sizes
        N = valuesN(nID);
        T = valuesT(tID);
        
        % Temporary arrays for results at each grid point (theta_{i1})    
        msePointStructsArray = cell(numTargetPoints, 1);   
        biasPointStructsArray = cell(numTargetPoints, 1);  
        varPointStructsArray = cell(numTargetPoints, 1);  
        weightsRegArray = cell(numTargetPoints, 1);  
        averageFirstWeightArray = cell(numTargetPoints, 1);
        averageMaxDiffArray = cell(numTargetPoints, 1);
        averageMassDiffArray = cell(numTargetPoints, 1);
        unitsUnresrtArray = cell(numTargetPoints, 1);
        
        % Loop through the grid of candidate theta_{i1} values
        for targetValueID = 1:numTargetPoints
            
            % Extract current target value
            targetValue = theta1Range(targetValueID);
            
            % Preallocate storage for replications
            errorsArrayTarget = cell(numReplicationsMSE, 1);
            estArrayTarget = errorsArrayTarget;
            targetParamsPoint  = errorsArrayTarget;
            
            % If weights are to be saved, create suitable arrays too
            if saveWeights 
                thetaPointArray = errorsArrayTarget;
                weightsArray = errorsArrayTarget;
                unitsUnrestrArray = errorsArrayTarget;
            end
            
            % Draw samples with current target value
            parfor replID=1:numReplicationsMSE 
                
                % --- Data Generation ---
                % Draw coefficients and error-term variances
                [thetaTrue, sigmaSq, thetaClassLabels] = ...
                        drawCoefficients(...
                            coefApproach, N, ...
                            varNoiseVar,  replID);
                 
                % Modify the first unit to have theta_{i1} = targetValue
                thetaTrue(1, 1) = targetValue;              
                
                % Simulate data
                [y, covars, u] = ...
                    drawData(thetaTrue, sigmaSq, varianceX, ...
                    T, replID);
                
                % Compute the true values of parameters for unit 1
                targetParamsPoint{replID} = computeAllParams(...
                    paramArray, thetaTrue(:, 1), ...
                    covars(:, 1, :), y(:, 1));
                
                % --- Individual Estimation ---
                % Estimate coefficients and variances using OLS
                [thetaHat, estVarianceArray] = OLS(y, covars);
                
                % --- Unit Averaging ---
                % Create optimal averaging schemes for the sample
                optimalSchemes = ...
                    createOptimalSchemes(...
                        thetaHat, thetaTrue, thetaClassLabels, ...
                        'firstOnly', ...
                        averagingIncludeBool);
                
                % Preprocess the methods array: expand generic optimal
                % position to use restrictions imposed by optimalSchemes  
                preparedMethodsArray = ...
                    addOptimalToMethodsStruct(methodsArray, optimalSchemes);
                
                % Apply averaging, compute errors, weights, unrestr. units
                [errorStruct, estStruct, weightStruct, unitsUnrestrStruct] = ...
                    sampleAveraging(y, covars,...
                    thetaHat,thetaTrue, estVarianceArray, ...
                    paramArray, preparedMethodsArray, ...
                    'firstOnly');
                
                % Save replication results
                errorsArrayTarget{replID} = errorStruct;
                estArrayTarget{replID} = estStruct;
                if saveWeights 
                    thetaPointArray{replID} = thetaTrue;
                    weightsArray{replID} = weightStruct; 
                    unitsUnrestrArray{replID} = unitsUnrestrStruct;
                end
                
                % Display progress information
                fprintf('N=%d, T=%d, Target: %d/%d, Replication %d \n', ...
                    N, T, targetValueID, numTargetPoints, replID);
            end
            
            % --- Process Results for Current Grid Point ---
            % MSE, bias, variance estimation
            [msePointStructsArray{targetValueID},...
               biasPointStructsArray{targetValueID},...
               varPointStructsArray{targetValueID}] = ...
                processOneGridValueMSE(...
                errorsArrayTarget, estArrayTarget, ...
                paramArray, targetParamsPoint);
            
            % Process the weights and unrestricted units if these are saved
            if saveWeights
                [weightsRegArray{targetValueID}, ...
                    averageFirstWeightArray{targetValueID}, ...
                    averageMaxDiffArray{targetValueID}, ...
                    averageMassDiffArray{targetValueID}, ...
                    unitsUnresrtArray{targetValueID}] = ...
                    processOneGridValueWeightsUnrestUnits(...
                        paramArray, thetaPointArray,...
                        weightsArray, unitsUnrestrArray, ...
                        theta1Range, 0.1);
            end
        end

        % --- Aggregate Results for All Grid Points ---
        % Process MSE, bias, and variance
        mseTablesNT{nID, tID} = ...
             combineAllPointsTable(msePointStructsArray, paramArray);
        biasTablesNT{nID, tID} = ...
             combineAllPointsTable(biasPointStructsArray, paramArray);      
        varTablesNT{nID, tID} = ...
             combineAllPointsTable(varPointStructsArray, paramArray);
         
        % Process weights and unrestricted units results if necessary
        if saveWeights
            weightsTablesNT{nID, tID} = ...
                combineAllPointsWeightsUnitsTables(...
                    weightsRegArray, paramArray);
            firstWeightNT{nID, tID} = ...
                combineAllPointsWeightsUnitsTables(...
                    averageFirstWeightArray, paramArray);
            unitsUnrestrNT{nID, tID} = ...
                combineAllPointsWeightsUnitsTables(...
                    unitsUnresrtArray, paramArray);
            maxDiffNT{nID, tID} = ...
                combineAllPointsWeightsUnitsTables(...
                    averageMaxDiffArray, paramArray);
            massDiffNT{nID, tID} = ...
                combineAllPointsWeightsUnitsTables(...
                    averageMassDiffArray, paramArray);
        end
    end
end
 
%% Save and Export Results
close all; % Close any open figures
fileSaveName = makeOutputFileName('MSE', simulationSetting, ...
                                  numReplicationsMSE, valuesN, valuesT);
save(fileSaveName); % Save simulation results