%% Simulation auxiliary computations
% Obtain



T=valuesT;


% Set range for parameter to be explored
% eta1range = 0:0.03:0.99;
eta1range = 0:0.03*sqrt(T):sqrt(T)*0.99;
% eta1range = 0:0.03*sqrt(T):sqrt(T)*0.5;

etaRangeLen = length(eta1range);


% Obtain lengths
maxN = max(valuesN);
numN = length(valuesN);
numParams = length(mu); % obtain number of parameters used
meanCoef = [lambdaMean; meanBeta];
% Draw coefficients and variances





%% Main Loop
% Loop through parameter vector, draw multiple samples for each
% value


% Create MSE arrays
uaMSEarrays

% Main loop: loop through the range of eta1
for i=1:etaRangeLen
    
    % Change coordinate and obtain deviations
    %    thetaLoop(1, 1) =    meanCoef(1)+lambda1range(i)/sqrt(T);
    
    % Recreate temporary variance vectors
    uaLoopArrays
    
    % Inner loop: drawing samples
    for j=1:numReplications % parfor this
        
        %%%%%%%%%%%%%%%%%%%%%%%
        %%% Data Generation %%%
        %%%%%%%%%%%%%%%%%%%%%%%
        
        if localAsy==1 % Local asymptotics
            [~, thetaSample, sigmaSq] = uaDrawCoefficients(maxN, T, ...
                meanCoef, varianceBeta, varNoiseVar,  j);
            
        else % Fixed parameter asymptotics, for use with large T
            [thetaSample,~ , sigmaSq] = uaDrawCoefficients(maxN, 1, ...
                meanCoef, varianceBeta, varNoiseVar,  j);
        end
        
        a =    meanCoef(1)+eta1range(i)/sqrt(T);
        thetaLoop = thetaSample;
        thetaLoop(1, 1)= a;
        etaTrueLoop = sqrt(T)*(thetaLoop - repmat(meanCoef, 1, maxN)); % obtain deviations from the mean
        
        % Draw data, estimate coefficients and variances
        [y, covars, u] = uaSimulateData(thetaLoop, sigmaSq, varianceX,T, j);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Individual Estimation %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Estimate coefficients and variances
        [thetaHat, estVarianceArray] = uaOLS(y, covars);
        
        
        
        % Estimate variance
        Vest = T*linearDynamicVarianceEstimator(y, covars, thetaHat);
        % Form eta estimators
        etaEst = sqrt(T)*(thetaHat - repmat(mean(thetaHat,2),1, maxN));
        
        % IMPLEMENTATION NOTE:
        % Vest holds asymptotic variances. Accordingly, we multiply etaEst
        % by sqrt(T). This approach explicitly agrees with the notation in
        % the paper. However, we can drop multiplication by T without
        % changing the result
        
        
        %%%%%%%%%%%%%%%%%%%%%%
        %%% Unit Averaging %%%
        %%%%%%%%%%%%%%%%%%%%%%
        
        % Loop over parameters of interest
        for parID=1:numParams
            
            % Individual estimator of parameters
            unitParamEsts = mu{parID}(thetaHat, covars, y);
            % True target parameter value
            targetParamValues = mu{parID}(thetaLoop(:, 1), covars(:, 1,:), y(:, 1));
            % Estimate gradient of mu
            gradientEst  = D{parID}(thetaHat, covars, y);
            
            
            % Fixed-N weights
            uaAveragingWeights(thetaHat, estVarianceArray, gradientEst, 1)
            
            for unit = 1:3
               1; 
            end
            
            % Fixed-N weights and large-N weights
            [errInLoopPlugInFixed(parID, j,:), ...
                errInLoopPlugInLargeRandom1(parID, j,:),...
                errInLoopPlugInLargeSmallBias1(parID, j,:), ...
                errInLoopPlugInLargeLargeBias1(parID, j, :)] = ...
                uaEstimationErrorInfeasibleWeights(etaEst, V, gradientEst, ...
                i0_1, valuesN, targetParamValues, unitParamEsts);
            [~, ...
                errInLoopPlugInLargeRandom2(parID, j,:),...
                errInLoopPlugInLargeSmallBias2(parID, j,:),...
                errInLoopPlugInLargeSmallBias2(parID, j, :)] = ...
                uaEstimationErrorInfeasibleWeights(etaEst, V, gradientEst, ...
                i0_2, valuesN, targetParamValues, unitParamEsts);
            
            
            % Infeasible weights: using the infeasible psi matrix
            % Observe that the function receives true etas, variance and
            % true population gradient
            % These infeasible weights can be used to check the sanity of
            % the approximation to the MSE derived to the paper. The
            % optimized MSE returned by theses weights should also be no
            % worse than that of the individual estimator
            [errInLoopInfeasiblePsiFixed(parID, j,:), ...
                errInLoopInfeasiblePsiLargeRandom1(parID, j,:),...
                errInLoopInfeasiblePsiLargeSmallBias1(parID, j,:), ...
                errInLoopInfeasiblePsiLargeLargeBias1(parID, j, :)] =...
                uaEstimationErrorInfeasibleWeights(etaTrueLoop, V, d0,...
                i0_1, valuesN, targetParamValues, unitParamEsts);
            
            [~, errInLoopInfeasiblePsiLargeRandom2(parID, j,:),...
                errInLoopInfeasiblePsiLargeSmallBias2(parID, j,:), errInLoopInfeasiblePsiLargeSmallBias2(parID, j, :)] =...
                uaEstimationErrorInfeasibleWeights(etaTrueLoop, V, d0, i0_2, valuesN, targetParamValues, unitParamEsts);
            
            % AIC, MMA weights
            [errInLoopAIC(parID, j,:), errInLoopMMA(parID, j,:)] =...
                uaEstimationErrorAICMMAWWeights(thetaHat, y, covars, valuesN, targetParamValues, unitParamEsts); %#ok<*SAGROW>
            
            % Individual
            errInLoopIndividual(parID, j) = unitParamEsts(1)-targetParamValues;
            
            % Mean group
            errInLoopMG(parID, j,:) = uaEstimationErrorMeanGroup(valuesN, targetParamValues,unitParamEsts);
            
            
        end
        [j, i]
    end
    
    uaMSEloop
end


fileSaveName =   "Outputs/"+num2str(numReplications)+ "etaRange"+...
    num2str(min(eta1range))+ "-" + num2str(ceil(max(eta1range)))+ "N"+...
    num2str(valuesN(end))+"T"+num2str(T)+".mat";

save(fileSaveName)
