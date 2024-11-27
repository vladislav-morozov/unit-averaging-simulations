% ===========================================================
% File: simulateCI.m
% Description: This script simulates coverage and length properties of
%              confidence intervals (CIs) derived from unit averaging
%              estimators in heterogeneous panel data models.
% ===========================================================
%
% Project Name: Unit Averaging for Heterogeneous Panels
% Authors: Christian Brownlees, Vladislav Morozov
%
% Model: Heterogeneous panel ARX(1) model
%   y_{it} = theta_{i1} y_{it-1} + theta_{i2} x_{it} + u_{it}.
%
% This file evaluates the performance of confidence intervals (CIs) in 
% terms of coverage and length. Specifically, it examines the properties 
% of CIs generated using bootstrapped estimates of unit averaging methods.
%
% Simulation Workflow:
% 1. For each combination of N and T (number of units, time periods):
%    a. Iterate over candidate values of theta_{i1} (grid points).
%      I. Iterate over samples for each candidate value:
%           i. Generate data.
%           ii. Bootstrap the dataset and compute the bootstrap CI.
%           iii. Check coverage and length of the bootstrap CI.
%      II. Aggregate results for the grid point (e.g., compute MSE).
%    b. Aggregate results over grid into a common format (see below).
% 2. Save and export results.
%
% Output format:
% - Results are saved in cell arrays indexed by combinations of (N, T)
% - Results for coverage and length are saved as  Tables with columns
%  indexing unit averaging approaches.
% ===========================================================

%% Simulation Parameters and Initialization

% Extract dimensions of the vectors defining the simulation
numT = length(valuesT); % Number of time periods
numN = length(valuesN); % Number of cross-sectional units
numParams = length(paramArray); % Number of parameters considered
numTargetPoints = length(theta1Range); % Grid points for target theta_{i1}
 
% Preallocate results arrays for CI coverage and length
coverageNT = cell(numN, numT);
lengthNT = cell(numN, numT);

%% Main simulation loop
 
% Loop over the number of time periods (T) and number of units (N)
for tID = 1:numT
    for nID=1:numN
        % Current sample sizes
        N = valuesN(nID);
        T = valuesT(tID);
        
        % Temporary arrays for results at each grid point (theta_{i1})
        coveragePointsArray = cell(numTargetPoints, 1);
        lengthPointsArray = coveragePointsArray;
        
        % Loop through the grid of candidate theta_{i1} values
        for targetValueID = 1:numTargetPoints
            
            % Extract target value for first unit
            targetValue = theta1Range(targetValueID);
            
            % Preallocate storage for each replication
            coverageArrayTarget = cell(numReplicationsCI, 1);
            lengthArrayTarget = coverageArrayTarget;
            
            % Draw samples with current target value
            for replID=1:numReplicationsCI
                
                % --- Data Generation ---
                % Draw coefficients and error-term variances
                [thetaTrue, sigmaSq, thetaLabels] = ...
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
                targetParams = computeAllParams(...
                    paramArray, thetaTrue(:, 1), ...
                    covars(:, 1, :), y(:, 1));
                
                % --- Bootstrapping ---
                % Preallocate array to hold boostrap estimates
                bsEstsArrayTarget = cell(numBootstrapSamples, 1);

                % Generate bootstrap samples
                for bsSampleID = 1:numBootstrapSamples % parfor this
                    % Create bootstrap sample with the stationary bootstrap
                    [yBs,covarsBs] = ...
                        stationaryBootstrap(y, covars, ceil(T^(1/3)), T);
                    
                    % Estimate parameters on bootstrap sample
                    [thetaHatBs, estVarianceArrayBs] = OLS(yBs, covarsBs);
                    
                    % Create optimal averaging strategies
                    optimalSchemes = ...
                        createOptimalSchemes(...
                        thetaHatBs, thetaTrue, thetaLabels, ...
                        'firstOnly', ...
                        averagingIncludeBool);
                    preparedMethodsArray = ...
                        addOptimalToMethodsStruct(methodsArray, optimalSchemes);
                    
                    % Apply unit averaging and save estimates
                    [~, estStructBs, weightStructBs] = ...
                        sampleAveraging(yBs, covarsBs,...
                        thetaHatBs,thetaTrue,...
                        estVarianceArrayBs, paramArray, preparedMethodsArray, ...
                        'firstOnly');
                    
                    % Compute bias-corrected bootstrap estimates
                    estDebiasedStructBs = ...
                        debiasedBootstrapEstimator(...
                            estStructBs, weightStructBs, ...
                            paramArray, ...
                            thetaHatBs, covarsBs, yBs);
                    
                    % Save estimates
                    bsEstsArrayTarget{bsSampleID} = estDebiasedStructBs;
                end
                
                % --- Process Results for Current Sample --- 
                % Process bootstrap estimates into a table
                bsEstTableStruct = ...
                    processOneSampleBootstrapEstimates(...
                    bsEstsArrayTarget, paramArray...
                    );
                
                % Compute coverage and length of CIs in this dataset
                [~, coverStruct, lengthStruct] = ...
                    computeIntervalCoverageLength(...
                    bsEstTableStruct, paramArray, targetParams, alphaCI);
                
                % Store coverage and length results
                coverageArrayTarget{replID} = coverStruct;
                lengthArrayTarget{replID} = lengthStruct;
                
                % Print iteration information
                fprintf('N=%d, T=%d, Target: %d/%d, Replication %d \n', ...
                    N, T, targetValueID, numTargetPoints, replID);
            end
            
            % --- Process Results for Current Grid Point ---
            % Compute coverage for this value
            [coveragePointsArray{targetValueID}, ...
                lengthPointsArray{targetValueID}] = ...
                processOneGridValueCoverageLength(coverageArrayTarget, ...
                lengthArrayTarget, paramArray); 
        end
        
        % --- Aggregate Results for All Grid Points ---
        % Combine coverages and lengths into a table
        coverageNT{nID, tID} = ...
             combineAllPointsTable(coveragePointsArray, paramArray);     
        lengthNT{nID, tID} = ...
             combineAllPointsTable(lengthPointsArray, paramArray);
    end
end

%% Save and Export Results
fileSaveName = makeOutputFileName('CI', simulationSetting, ...
                                  numReplicationsCI, valuesN, valuesT);
save(fileSaveName); % Save simulation results