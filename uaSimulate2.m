



% Set range for parameter to be explored
% eta1range = 0:0.03:0.99;
% eta1range = 0:0.03*sqrt(T):sqrt(T)*0.99;
eta1range = 0:0.03*sqrt(T):sqrt(T)*0.5;

endOfRange = length(eta1range);


% Obtain lengths
N = max(Nvalues);
Nlen = length(Nvalues);
numPar = length(mu); % obtain number of parameters used
meanCoef = [lambdaMean; meanBeta];
% Draw coefficients and variances



 

% Main loop: loop through parameter vector, draw multiple samples for each
% value
 

% Create MSE arrays
uaMSEarrays
for i=1:endOfRange
    
   % Change coordinate and obtain deviations
%    thetaLoop(1, 1) =    meanCoef(1)+lambda1range(i)/sqrt(T);

   % Recreate temporary variance vectors
   uaLoopArrays
   
   % Inner loop: drawing samples
    parfor j=1:numReplications
        
        if localAsy==1 % Local asymptotics
            [~, thetaSample, sigmaSq] = uaDrawCoefficients(N, T, ...
                meanCoef, varianceBeta, varNoiseVar,  j);
            
        else % Fixed parameter asymptotics, for use with large T
            [thetaSample,~ , sigmaSq] = uaDrawCoefficients(N, 1, ...
                meanCoef, varianceBeta, varNoiseVar,  j);
        end

        a =    meanCoef(1)+eta1range(i)/sqrt(T);
        thetaLoop = thetaSample;
        thetaLoop(1, 1)= a;
        etaTrueLoop = sqrt(T)*(thetaLoop - repmat(meanCoef, 1, N)); % obtain deviations from the mean
        % Obtain true values of variance
        V = uaTrueAsymptoticVariance(thetaSample(1,:), thetaSample(2,:),sigmaSq, N,1);    % compute population-fixed variances
        
        % Draw data, estimate coefficients and variances
        [y, x, u] = uaSimulateData(thetaLoop, sigmaSq, varianceX,T, j); % draw data
        thetaHat = linearStaticEstimators(y, x);% estimate coefficients
        Vest = T*linearDynamicVarianceEstimator(y, x, thetaHat); % estimate variance
        etaEst = sqrt(T)*(thetaHat - repmat(mean(thetaHat,2),1, N));
        

        % Loop over parameters
        for par=1:numPar
            % Gradient, slice out of the supplied array of funcitons
            d0  = D{par}(thetaLoop, x, y); % Gradient describes which parameter is estimated
            d1  = D{par}(thetaHat, x, y); % No estimation needed in this case
            
            % Target value
            targetValue = mu{par}(thetaLoop(:, 1), x(:, 1,:), y(:, 1)); % set target value
            unitEst = mu{par}(thetaHat, x, y);
            
            % Infeasible weights: using the infeasible psi matrix
            [errInLoopInfeasiblePsiFixed(par, j,:), errInLoopInfeasiblePsiLargeRandom1(par, j,:),...
                errInLoopInfeasiblePsiLargeSmallBias1(par, j,:), errInLoopInfeasiblePsiLargeLargeBias1(par, j, :)] =...
                uaEstimationErrorInfeasibleWeights(etaTrueLoop, V, d0, i0_1, Nvalues, targetValue, unitEst); %#ok<*SAGROW>
            [~, errInLoopInfeasiblePsiLargeRandom2(par, j,:),...
                errInLoopInfeasiblePsiLargeSmallBias2(par, j,:), errInLoopInfeasiblePsiLargeSmallBias2(par, j, :)] =...
                uaEstimationErrorInfeasibleWeights(etaTrueLoop, V, d0, i0_2, Nvalues, targetValue, unitEst);

            % Feasible plug-in weights
            [errInLoopPlugInFixed(par, j,:), errInLoopPlugInLargeRandom1(par, j,:),...
                errInLoopPlugInLargeSmallBias1(par, j,:), errInLoopPlugInLargeLargeBias1(par, j, :)] =...
                uaEstimationErrorInfeasibleWeights(etaEst, V, d1, i0_1, Nvalues, targetValue, unitEst); %#ok<*SAGROW>
            [~, errInLoopPlugInLargeRandom2(par, j,:),...
                errInLoopPlugInLargeSmallBias2(par, j,:), errInLoopPlugInLargeSmallBias2(par, j, :)] =...
                uaEstimationErrorInfeasibleWeights(etaEst, V, d1, i0_2, Nvalues, targetValue, unitEst);

            % AIC, MMA weights
            [errInLoopAIC(par, j,:), errInLoopMMA(par, j,:)] =...
                 uaEstimationErrorAICMMAWWeights(thetaHat, y, x, Nvalues, targetValue, unitEst); %#ok<*SAGROW>
            
            % Individual 
            errInLoopIndividual(par, j) = unitEst(1)-targetValue;    

            % Mean group
            errInLoopMG(par, j,:) = uaEstimationErrorMeanGroup(Nvalues, targetValue,unitEst);
 
          
        end
        [j, i]
    end
    uaMSEloop
end


fileSaveName =   "Outputs/"+num2str(numReplications)+ "etaRange"+...
        num2str(min(eta1range))+ "-" + num2str(ceil(max(eta1range)))+ "N"+...
        num2str(Nvalues(end))+"T"+num2str(T)+".mat";
    
save(fileSaveName)
