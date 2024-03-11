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
    parfor j=1:numReplications
        
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
        % Obtain true values of variance
        V = uaTrueAsymptoticVariance(thetaSample(1,:), thetaSample(2,:),sigmaSq, maxN,1);    % compute population-fixed variances
        
        % Draw data, estimate coefficients and variances
        [y, x, u] = uaSimulateData(thetaLoop, sigmaSq, varianceX,T, j); % draw data
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Individual Estimation %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Estimate coefficients
        thetaHat = linearStaticEstimators(y, x);%
        % Estimate variance
        Vest = T*linearDynamicVarianceEstimator(y, x, thetaHat);
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
            
            % Compute the gradient of mu
            % True gradient (for infeasible weights)
            d0  = D{parID}(thetaLoop, x, y);
            % Estimated gradient (for feasible weights)
            d1  = D{parID}(thetaHat, x, y);
            
            % True target parameter value
            targetValue = mu{parID}(thetaLoop(:, 1), x(:, 1,:), y(:, 1));
            
            % Individual estimates of parameter of interest
            unitEst = mu{parID}(thetaHat, x, y);
            
            % Fixed-N weights and large-N weights
            [errInLoopPlugInFixed(parID, j,:), ...
                errInLoopPlugInLargeRandom1(parID, j,:),...
                errInLoopPlugInLargeSmallBias1(parID, j,:), ...
                errInLoopPlugInLargeLargeBias1(parID, j, :)] = ...
                    uaEstimationErrorInfeasibleWeights(etaEst, V, d1, ...
                     i0_1, valuesN, targetValue, unitEst); 
            [~, ...
                errInLoopPlugInLargeRandom2(parID, j,:),...
                errInLoopPlugInLargeSmallBias2(parID, j,:),...
                errInLoopPlugInLargeSmallBias2(parID, j, :)] = ...
                    uaEstimationErrorInfeasibleWeights(etaEst, V, d1, ...
                    i0_2, valuesN, targetValue, unitEst);
            
            
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
                        i0_1, valuesN, targetValue, unitEst);
            
            [~, errInLoopInfeasiblePsiLargeRandom2(parID, j,:),...
                errInLoopInfeasiblePsiLargeSmallBias2(parID, j,:), errInLoopInfeasiblePsiLargeSmallBias2(parID, j, :)] =...
                uaEstimationErrorInfeasibleWeights(etaTrueLoop, V, d0, i0_2, valuesN, targetValue, unitEst);
            
            % AIC, MMA weights
            [errInLoopAIC(parID, j,:), errInLoopMMA(parID, j,:)] =...
                uaEstimationErrorAICMMAWWeights(thetaHat, y, x, valuesN, targetValue, unitEst); %#ok<*SAGROW>
            
            % Individual
            errInLoopIndividual(parID, j) = unitEst(1)-targetValue;
            
            % Mean group
            errInLoopMG(parID, j,:) = uaEstimationErrorMeanGroup(valuesN, targetValue,unitEst);
            
            
        end
        [j, i]
    end
    
    uaMSEloop
end


fileSaveName =   "Outputs/"+num2str(numReplications)+ "etaRange"+...
    num2str(min(eta1range))+ "-" + num2str(ceil(max(eta1range)))+ "N"+...
    num2str(valuesN(end))+"T"+num2str(T)+".mat";

save(fileSaveName)
