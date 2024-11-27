%% Simulation auxiliary computations

% Obtain lengths
numT = length(valuesT);
numN = length(valuesN);
numParams = length(paramArray); % obtain number of parameters used


numTargetPoints = length(theta1Range);

% Implementation
% Loop through (samples), (N, T) pairs
% Save three cell arrays with those indices
% 1. trueTheta array -- hold true values
% 2. weight array -- deeper level ->.parID -> .schemeID -> N_kx1 vector of
% errors. jth entry corresponds to par(j)th unit being the target. Matches
% up with the rows of the trueThetas.
% 2. error array -- deeper level ->.parID -> .schemeID -> .weights and
% .unrestricted. Both house N_kxN_k matrices, jth rows correspond to the
% jth unit being the targe
% The inloop functions return those deeper levels


%% Main Loop
% Loop through parameter vector, draw multiple samples for each
% value

coverageNT = cell(numN, numT);
lengthNT = cell(numN, numT);

% Loop through the different values of N
for tID = 1:numT
    for nID=1:numN
        % Current sample sizes
        N = valuesN(nID);
        T = valuesT(tID);
        
        % Create array to hold coverage
        coveragePointsArray = cell(numTargetPoints, 1);
        lengthPointsArray = coveragePointsArray;
        
        % Iterate through the target values
        for targetValueID = 1:numTargetPoints
            
            % Extract target value for first unit
            targetValue = theta1Range(targetValueID);
            
            % Create arrays to fill out in the loop
            coverageArrayTarget = cell(numReplicationsCI, 1);
            lengthArrayTarget = coverageArrayTarget;
            
            % Draw samples with current target value
            for replID=1:numReplicationsCI
                
                %%%%%%%%%%%%%%%%%%%%%%%
                %%% Data Generation %%%
                %%%%%%%%%%%%%%%%%%%%%%%
                
                % Draw sample
                [thetaTrue, sigmaSq, thetaLabels] = ...
                    drawCoefficients(...
                    coefApproach, N, ...
                    varNoiseVar,  replID);
                
                % Modify the first unit in line with the targetValue
                thetaTrue(1, 1) = targetValue;
                
                % Draw data, estimate coefficients and variances
                [y, covars, u] = ...
                    drawData(thetaTrue, sigmaSq, varianceX, ...
                    T, replID);
                
                % Compute the true values of parameters for unit 1
                targetParams = computeAllParams(...
                    paramArray, thetaTrue(:, 1), ...
                    covars(:, 1, :), y(:, 1));
                
                %%%%%%%%%%%%%%%%%%%%%
                %%% Bootstrapping %%%
                %%%%%%%%%%%%%%%%%%%%%
                
                % Create array to hold boostrap estimates
                bsEstsArrayTarget = cell(numBootstrapSamples, 1);
                % Draw bootstrap samples
                for bsSampleID = 1:numBootstrapSamples % parfor this
                    % Draw boostrap sample
                    [yBs,covarsBs] = ...
                        stationaryBootstrap(y, covars, ceil(T^(1/3)), T);
                    
                    % Reestimate
                    [thetaHatBs, estVarianceArrayBs] = OLS(yBs, covarsBs);
                    
                    % Create optimal strategies for this subsample
                    optimalSchemes = ...
                        createOptimalSchemes(...
                        thetaHatBs, thetaTrue, thetaLabels, ...
                        'firstOnly', ...
                        averagingIncludeBool);
                    preparedMethodsArray = ...
                        addOptimalToMethodsStruct(methodsArray, optimalSchemes);
                    
                    % Apply averaging and save estimates
                    [~, estStructBs, weightStructBs] = ...
                        sampleAveraging(yBs, covarsBs,...
                        thetaHatBs,thetaTrue,...
                        estVarianceArrayBs, paramArray, preparedMethodsArray, ...
                        'firstOnly');
                    
                    % Compute the bias-corrected bootstrap values
                    estDebiasedStructBs = ...
                        debiasedBootstrapEstimator(...
                            estStructBs, weightStructBs, ...
                            paramArray, ...
                            thetaHatBs, covarsBs, yBs);
                    
                    % Save estimates
                    bsEstsArrayTarget{bsSampleID} = estDebiasedStructBs;
                    
                end
                
                bsEstTableStruct = ...
                    processOneSampleBootstrapEstimates(...
                    bsEstsArrayTarget, paramArray...
                    );
                
                % Obtain bootstrap coverage and length in this iteration
                [~, coverStruct, lengthStruct] = ...
                    computeIntervalCoverageLength(...
                    bsEstTableStruct, paramArray, targetParams, alphaCI);
                
                % Save the coverage and length structs
                coverageArrayTarget{replID} = coverStruct;
                lengthArrayTarget{replID} = lengthStruct;
                
                % Print iteration information
                fprintf('N=%d, T=%d, Target: %d/%d, Replication %d \n', ...
                    N, T, targetValueID, numTargetPoints, replID);
            end
            
            % Compute coverage for this value
            [coveragePointsArray{targetValueID}, ...
                lengthPointsArray{targetValueID}] = ...
                processOneGridValueCoverageLength(coverageArrayTarget, ...
                lengthArrayTarget, paramArray); 
        end
        
        % Combine coverages and lengths into a table
        coverageNT{nID, tID} = ...
             combineAllPointsTable(coveragePointsArray, paramArray);     
        lengthNT{nID, tID} = ...
             combineAllPointsTable(lengthPointsArray, paramArray);
    end
end

%% Save results

% Generate file name
fileSaveName = makeOutputFileName('CI', simulationSetting, numReplicationsCI, valuesN, valuesT);

% Save file
save(fileSaveName)