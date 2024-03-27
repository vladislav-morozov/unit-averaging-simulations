%

% for (N, T)
%   for candidate values of (theta_{i1}, theta_{i2})
%       for sample
%           draw thetas and set first unit to candidate value
%           simulate data
%           estimate parameters
%           conduct averaging on first unit
%           record errors of approaches
%       estimate MSE value for candidate value
%

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

mseTablesN = cell(numN, numT);
biasTablesN = cell(numN, numT);
varTablesN = cell(numN, numT);
weightsTablesNT = cell(numN, numT);
firstWeightNT = cell(numN, numT);
 
% Loop through the different values of N
for tID = 1:numT
    for nID=1:numN
        % Current sample sizes
        N = valuesN(nID);
        T = valuesT(tID);
        
        msePointStructsArray = cell(numTargetPoints, 1);   
        biasPointStructsArray = cell(numTargetPoints, 1);  
        varPointStructsArray = cell(numTargetPoints, 1);  
        weightsRegArray = cell(numTargetPoints, 1);  
        averageFirstWeightArray = cell(numTargetPoints, 1);  
        % Iterate through the target values
        for targetValueID = 1:numTargetPoints
            
            targetValue = theta1Range(targetValueID);
            
            % Create arrays to fill out in the loop
            errorsArrayTarget = cell(numReplications, 1);
            estArrayTarget = errorsArrayTarget;
            targetParamsPoint  = errorsArrayTarget;
            
            if saveWeights
                % Save sample and corresponding weights
                thetaPointArray = errorsArrayTarget;
                weightsArray = errorsArrayTarget;
            end
            
            % Draw samples with current target value
            parfor replID=1:numReplications % parfor this
                
                %%%%%%%%%%%%%%%%%%%%%%%
                %%% Data Generation %%%
                %%%%%%%%%%%%%%%%%%%%%%%
                
 
                [thetaTrue, sigmaSq, thetaLabels] = ...
                        uaDrawCoefficients(...
                            coefApproach, N, ...
                            varNoiseVar,  replID);
                 
                
                % Modify the first unit in line with the targetValue
                thetaTrue(1, 1) = targetValue;              
                
                % Draw data, estimate coefficients and variances
                [y, covars, u] = ...
                    uaSimulateData(thetaTrue, sigmaSq, varianceX, ...
                    T, replID);
                
                % Compute the true values of parameters for unit 1
                targetParamsPoint{replID} = uaComputeAllParams(...
                    paramArray, thetaTrue(:, 1), ...
                    covars(:, 1, :), y(:, 1));
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%% Individual Estimation %%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                % Estimate coefficients and variances
                [thetaHat, estVarianceArray] = uaOLS(y, covars);
                
                
                %%%%%%%%%%%%%%%%%%%%%%
                %%% Unit Averaging %%%
                %%%%%%%%%%%%%%%%%%%%%%
                
                % Create optimal strategies for this subsamples
                optimalSchemes = ...
                    uaOptimalSchemes(...
                        thetaHat, thetaTrue, thetaLabels, ...
                        'firstOnly', ...
                        averagingIncludeBool);
                
                % Preprocess the methods array: expand generic optimal position to
                % use the restrictions imposed by optimalSchemes
                preparedMethodsArray = ...
                    uaAddOptimalToMethodsStruct(methodsArray, optimalSchemes);
                
                % Apply averaging
                [errorStruct, estStruct, weightStruct, unitsUnrestrStruct] = ...
                    uaSampleAveraging(y, covars,...
                    thetaHat,thetaTrue,...
                    estVarianceArray, paramArray, preparedMethodsArray, ...
                    'firstOnly');
                
                % Save weights and errors
                errorsArrayTarget{replID} = errorStruct;
                estArrayTarget{replID} = estStruct;
                if saveWeights
                    % Save sample and corresponding weights
                    thetaPointArray{replID} = thetaTrue;
                    weightsArray{replID} = weightStruct;
                end
                if saveUnrestricted
                    unitsUnrestrArray = unitsUnrestrStruct;
                end
                
                
                disp(['N=', num2str(N), ...
                      ', T=', num2str(T), ...
                      ', Target: ', num2str(targetValueID),...
                      '/', num2str(numTargetPoints),...
                      ', Replication ', num2str(replID)])
            end
            
            % Compute MSE for current target value
            [msePointStructsArray{targetValueID},...
               biasPointStructsArray{targetValueID},...
               varPointStructsArray{targetValueID}] = ...
                uaProcessOneValueError(...
                errorsArrayTarget, estArrayTarget, ...
                paramArray, targetParamsPoint);
            
            % Process the weights if necessary
            if saveWeights
                [weightsRegArray{targetValueID}, ...
                    averageFirstWeightArray{targetValueID}] = ...
                    uaAverageWeight(paramArray, thetaPointArray,...
                    weightsArray, ...
                    theta1Range, 0.1);
            end
        end
        % Glue together the msePointStructsArray into a struct of tables
        [mseTablesN{nID, tID}, ...
            biasTablesN{nID, tID}, ...
            varTablesN{nID, tID}] = ...
            uaProcessMSE(msePointStructsArray, biasPointStructsArray,...
            varPointStructsArray, paramArray);
        
        if saveWeights
            % Glue together weight arrays
            [weightsTablesNT{nID, tID}, firstWeightNT{nID, tID}] = ...
                uaProcessWeights(paramArray, ...
                weightsRegArray, ...
                averageFirstWeightArray...
                );
        end
    end
end
%%
% Extract range of methods for plotting
optimalSchemes = ...
    uaOptimalSchemes(randn(2, N), randn(2, N), randn(N, 1), ...
    'firstOnly', averagingIncludeBool);
allMethodsArray = ...
    uaAddOptimalToMethodsStruct(methodsArray, optimalSchemes);


%% Save results

titleN = "";
for nID=1:numN
    titleN = titleN + "-" + valuesN(nID) ;
end
titleT = "";
for tID=1:numT
    titleT = titleT + "-" + valuesT(tID);
end

fileSaveName =   "Outputs/"+...
    simulationSetting + "/" + ...
    "Replication-" + num2str(numReplications)+ ...
    "N" + titleN + ...
    "T" + titleT +...
    "weights" + num2str(saveWeights) + ...
    "unrestricted" + num2str(saveUnrestricted) + ...
    ".mat";

% fileSaveName =   "Outputs/"+num2str(numReplications)+ "etaRange"+...
%     num2str(min(eta1range))+ "-" + num2str(ceil(max(eta1range)))+ "N"+...
%     num2str(valuesN(end))+"T"+num2str(T)+".mat";

save(fileSaveName)
