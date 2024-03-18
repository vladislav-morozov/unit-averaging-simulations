%% Simulation auxiliary computations
% Obtain



T=valuesT;


% Obtain lengths
maxN = max(valuesN);
numN = length(valuesN);
numParams = length(paramArray); % obtain number of parameters used
meanCoef = [lambdaMean; meanBeta];


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

trueThetaArray = cell(numN, numReplications);
errorsArray = cell(numN, numReplications);
if saveWeights
    weightsArray = cell(numN, numReplications);
end
if saveUnrestricted
    unitsUnrestrArray = cell(numN, numReplications);
end

 
% Draw samples
parfor replID=1:numReplications % parfor this
    
    %%%%%%%%%%%%%%%%%%%%%%%
    %%% Data Generation %%%
    %%%%%%%%%%%%%%%%%%%%%%%
    
    localAsy = 0;
    if localAsy==1 % Local asymptotics
        [~, thetaTrue, sigmaSq] = uaDrawCoefficients(maxN, T, ...
            meanCoef, varianceBeta, varNoiseVar,  replID);
        
    else % Fixed parameter asymptotics, for use with large T
        [thetaTrue,~ , sigmaSq] = uaDrawCoefficients(maxN, 1, ...
            meanCoef, varianceBeta, varNoiseVar,  replID);
    end
    
    
    % Draw data, estimate coefficients and variances
    [y, covars, u] = uaSimulateData(thetaTrue, sigmaSq, varianceX,T, replID);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Individual Estimation %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Estimate coefficients and variances
    [thetaHat, estVarianceArray] = uaOLS(y, covars);
    
    
    %%%%%%%%%%%%%%%%%%%%%%
    %%% Unit Averaging %%%
    %%%%%%%%%%%%%%%%%%%%%%
    
    % Loop through the different values of N
    for nID=1:numN
        
        % Extract value of N
        N = valuesN(nID);
        
        % Draw N random units
        subsample = datasample(1:maxN, N);
        % Extract cross-sectional subsample
        thetaHatSub = thetaHat(:, subsample);
        thetaTrueSub = thetaTrue(:, subsample);
        estVarianceArraySub = estVarianceArray(subsample);
        ySub = y(:, subsample);
        covarsSub = covars(:, subsample, :);
        
        
        % Save thetas
        trueThetaArray{nID, replID} = thetaTrueSub';
        
        % Create optimal strategies for this subsamples
        optimalSchemes = uaOptimalSchemes(thetaHatSub);
        
        % Preprocess the methods array: expand generic optimal position to
        % use the restrictions imposed by optimalSchemes
        preparedMethodsArray = ...
            uaAddOptimalToMethodsStruct(methodsArray, optimalSchemes);
        
        % Apply averaging
        [errorStruct, weightStruct, unitsUnrestrStruct] = ...
            uaSampleAveraging(ySub, covarsSub,...
            thetaHatSub,thetaTrueSub,...
            estVarianceArraySub, paramArray, preparedMethodsArray);
        
        % Save weights and errors
        errorsArray{nID, replID} = errorStruct;
        if saveWeights
            weightsArray{nID, replID} = weightStruct;
        end
        if saveUnrestricted
            unitsUnrestrArray{nID, replID} = unitsUnrestrStruct;
        end
        
        
        
    end
    disp(['Replication ', num2str(replID)])
end

%%
titleN = "";
for nID=1:numN
    titleN = titleN + "-" + valuesN(nID) ;
end

fileSaveName =   "Outputs/"+...
    "design-1" + ...
    "Replication-" + num2str(numReplications)+ ...
    "N" + titleN + ...
    "T"+num2str(T)+...
    "weights" + num2str(saveWeights) + ...
    "unrestricted" + num2str(saveUnrestricted) + ...
    ".mat";

% fileSaveName =   "Outputs/"+num2str(numReplications)+ "etaRange"+...
%     num2str(min(eta1range))+ "-" + num2str(ceil(max(eta1range)))+ "N"+...
%     num2str(valuesN(end))+"T"+num2str(T)+".mat";

save(fileSaveName)
