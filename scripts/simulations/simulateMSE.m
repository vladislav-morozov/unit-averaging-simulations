% ===========================================================
% File: simulateMSE.m
% Description: This script runs the main simulation for evaluating the MSE
%   of unit averaging approaches

% ===========================================================
%
% Project Name: Unit Averaging for Heterogeneous Panels
% Authors: Christian Brownlees, Vladislav Morozov
%
% Model: Heterogeneous panel ARX(1) model
%   y_{it} = theta_{i1} y_{it-1} + theta_{i2} x_{it} + u_{it}
%
%
% This file simulates the MSE, weight, and unrestricted properties of the
% various unit averaging approaches. 
% Pseudocode for the simulation loop:
%
% for (N, T)
%   for candidate values of (theta_{i1}, theta_{i2})
%       for sample
%           draw thetas and set first unit to candidate value
%           simulate data
%           estimate parameters
%           conduct averaging on first unit
%           record errors of approaches
%       estimate MSE, bias, and variance for candidate value
%       if necessary, estimate weight properties
%   process results for all grid points 
% 
% Notes regarding format of results
% 1. Results are saved in cell arrays indexed by combinations of (N, T)
% 2. Results come in two formats
%     - Tables with columns indexing unit averaging approaches. Applies to 
%       MSE, bias, variance, differences in masses between unrestricted and
%       restricted units, average own weight
%     - Structs indexed by averaging approaches 

% ===========================================================
%% Extract dimensions of vectors that define the simulation

% Numbers of dimensions of samples drawn
numT = length(valuesT);
numN = length(valuesN);

% Number of target parameters
numParams = length(paramArray); 
 
% Number of points in the grid for \theta_{i1}
numTargetPoints = length(theta1Range);

%% Allocate arrays that will hold the results

% MSE
mseTablesNT = cell(numN, numT);
% Bias
biasTablesNT = cell(numN, numT);
% Variance
varTablesNT = cell(numN, numT);
% All weights
weightsTablesNT = cell(numN, numT);
% Weights assigned to first unit
firstWeightNT = cell(numN, numT);
% Maximum difference in weights
maxDiffNT = cell(numN, numT);
% Difference in mass between restricted and unrestricted sets
massDiffNT = cell(numN, numT); 
% Unrestricted units
unitsUnrestrNT = cell(numN, numT); 

%% Main Loop
 
for tID = 1:numT
    for nID=1:numN
        % Current sample sizes
        N = valuesN(nID);
        T = valuesT(tID);
        
        % Create temporary arrays to hold results for each in the grid for
        % theta_{i1}    
        msePointStructsArray = cell(numTargetPoints, 1);   
        biasPointStructsArray = cell(numTargetPoints, 1);  
        varPointStructsArray = cell(numTargetPoints, 1);  
        weightsRegArray = cell(numTargetPoints, 1);  
        averageFirstWeightArray = cell(numTargetPoints, 1);
        averageMaxDiffArray = cell(numTargetPoints, 1);
        averageMassDiffArray = cell(numTargetPoints, 1);
        unitsUnresrtArray = cell(numTargetPoints, 1);
        
        % Iterate through the target values
        for targetValueID = 1:numTargetPoints
            
            % Extract current target value
            targetValue = theta1Range(targetValueID);
            
            % Create arrays to fill out in the loop
            errorsArrayTarget = cell(numReplicationsMSE, 1);
            estArrayTarget = errorsArrayTarget;
            targetParamsPoint  = errorsArrayTarget;
            
            % If weights are to be saved, create suitable arrays to
            if saveWeights 
                thetaPointArray = errorsArrayTarget;
                weightsArray = errorsArrayTarget;
                unitsUnrestrArray = errorsArrayTarget;
            end
            
            % Draw samples with current target value
            for replID=1:numReplicationsMSE 
                
                %%%%%%%%%%%%%%%%%%%%%%%
                %%% Data Generation %%%
                %%%%%%%%%%%%%%%%%%%%%%%
                
                % Draw coefficients
                [thetaTrue, sigmaSq, thetaClassLabels] = ...
                        drawCoefficients(...
                            coefApproach, N, ...
                            varNoiseVar,  replID);
                 
                % Modify the first unit to have theta_{i1} = targetValue
                thetaTrue(1, 1) = targetValue;              
                
                % Draw data, estimate coefficients and variances
                [y, covars, u] = ...
                    drawData(thetaTrue, sigmaSq, varianceX, ...
                    T, replID);
                
                % Compute the true values of parameters for unit 1
                targetParamsPoint{replID} = computeAllParams(...
                    paramArray, thetaTrue(:, 1), ...
                    covars(:, 1, :), y(:, 1));
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%% Individual Estimation %%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                % Estimate coefficients and variances
                [thetaHat, estVarianceArray] = OLS(y, covars);
                
                %%%%%%%%%%%%%%%%%%%%%%
                %%% Unit Averaging %%%
                %%%%%%%%%%%%%%%%%%%%%%
                
                % Create optimal unit averaging schemes for this subsample
                optimalSchemes = ...
                    createOptimalSchemes(...
                        thetaHat, thetaTrue, thetaClassLabels, ...
                        'firstOnly', ...
                        averagingIncludeBool);
                
                % Preprocess the methods array: expand generic optimal position to
                % use the restrictions imposed by optimalSchemes
                preparedMethodsArray = ...
                    addOptimalToMethodsStruct(methodsArray, optimalSchemes);
                
                % Apply averaging
                [errorStruct, estStruct, weightStruct, unitsUnrestrStruct] = ...
                    sampleAveraging(y, covars,...
                    thetaHat,thetaTrue, estVarianceArray, ...
                    paramArray, preparedMethodsArray, ...
                    'firstOnly');
                
                % Save weights and errors
                errorsArrayTarget{replID} = errorStruct;
                estArrayTarget{replID} = estStruct;
                if saveWeights
                    % Save sample and corresponding weights
                    thetaPointArray{replID} = thetaTrue;
                    weightsArray{replID} = weightStruct;
                    % Save unrestricted units
                    unitsUnrestrArray{replID} = unitsUnrestrStruct;
                end
                
                % Print iteration information
                fprintf('N=%d, T=%d, Target: %d/%d, Replication %d \n', ...
                    N, T, targetValueID, numTargetPoints, replID);
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Processing for current grid point %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Compute MSE for current target value
            [msePointStructsArray{targetValueID},...
               biasPointStructsArray{targetValueID},...
               varPointStructsArray{targetValueID}] = ...
                processOneGridValueMSE(...
                errorsArrayTarget, estArrayTarget, ...
                paramArray, targetParamsPoint);
            
            % Process the weights if these are saved
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

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Processing for the whole DGP setting %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
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
                combineAllPointsWeightsUnitsTables(weightsRegArray, paramArray);
            firstWeightNT{nID, tID} = ...
                combineAllPointsWeightsUnitsTables(averageFirstWeightArray, paramArray);
            unitsUnrestrNT{nID, tID} = ...
                combineAllPointsWeightsUnitsTables(unitsUnresrtArray, paramArray);
            maxDiffNT{nID, tID} = ...
                combineAllPointsWeightsUnitsTables(averageMaxDiffArray, paramArray);
            massDiffNT{nID, tID} = ...
                combineAllPointsWeightsUnitsTables(averageMassDiffArray, paramArray);
        end
    end
end
 

%% Save results

% Close any open figures before exporting results
close all

% Create file name based on simulation parameters
fileSaveName = makeOutputFileName('MSE', simulationSetting, numReplicationsMSE, valuesN, valuesT);

% Export simulation results
save(fileSaveName)