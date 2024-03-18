%% Subfile for replicating figure for one-step forecasts
% Definiitions and declarations must carry over form main file

% mu is the one-step ahead forecasts for the dynamic file

% Purpose of this separate file is that the population matrix is
% constructed in the internal loop, since the Jacobian depends on a
% particular data set through y_{1,T}



% Recreate objects for later use
 % Create equal weights
 equalWeights1 = zeros(N,1);
 equalWeights1(1:compValue1) = ones(compValue1,1)/compValue1;
 equalWeights2 = zeros(N,1);
 equalWeights2(1:compValue2) = ones(compValue2,1)/compValue2;
 equalWeightsAll = ones(N,1)/(N);

 
 % Create vectors to house variances
    % Individual
 varianceIndividual = zeros(endOfRange, 1);
    % Plug-in optimized
 varianceSample1 = zeros(endOfRange, 1);
 varianceSample2 = zeros(endOfRange, 1);
 varianceSampleAll = zeros(endOfRange, 1);
    % Population optimal
 varianceP1 = zeros(endOfRange, 1);
 varianceP2 = zeros(endOfRange, 1);
 variancePAll = zeros(endOfRange, 1);
    % Stein-type
 varianceStein1 = zeros(endOfRange, 1);
 varianceStein2 = zeros(endOfRange, 1);
 varianceStein = zeros(endOfRange, 1);
    % Equal
 varianceEqual1 = zeros(endOfRange, 1);
 varianceEqual2 = zeros(endOfRange, 1);
 varianceEqualAll = zeros(endOfRange, 1);
    % AIC/BIC
 varianceAIC1 = zeros(endOfRange, 1);
 varianceAIC2 = zeros(endOfRange, 1);
 varianceAIC = zeros(endOfRange, 1);
    % MMA
 varianceMMA1 = zeros(endOfRange, 1);
 varianceMMA2 = zeros(endOfRange, 1);
 varianceMMA = zeros(endOfRange, 1);

[betaOriginal,~ , sigmaSq] = linearDynamicDrawCoefficients(N, T, design,... 
            meanBeta, varianceBeta, seedCoefficients);
% Main loop: loop through parameter vector, draw multiple samples
for j=1:endOfRange
   
   betaLoop =    dynamicChangeCoordinate1(lambdaRange(j), betaOriginal); 
   etaTrue = betaLoop - repmat(meanCoef, 1, N); % obtain deviations from the mean

   % Recreate variance vectors
   varianceInLoopSample1 = zeros(numReplications, 1);
   varianceInLoopSample2 = zeros(numReplications, 1);
   varianceInLoopSampleAll = zeros(numReplications, 1);
   varianceInLoopStein1 = zeros(numReplications, 1);
   varianceInLoopAIC1 = zeros(numReplications, 1);
   varianceInLoopAIC2 = zeros(numReplications, 1);
   varianceInLoopAIC = zeros(numReplications, 1);
   varianceInLoopMMA1 = zeros(numReplications, 1);
   varianceInLoopMMA2 = zeros(numReplications, 1);
   varianceInLoopMMA = zeros(numReplications, 1);
   varianceInLoopOptimalP1 = zeros(numReplications, 1);
   varianceInLoopOptimalP2 = zeros(numReplications, 1);
   varianceInLoopOptimalPAll = zeros(numReplications, 1);
   varianceInLoopEqual1 = zeros(numReplications, 1);
   varianceInLoopEqual2 = zeros(numReplications, 1);
   varianceInLoopEqualAll = zeros(numReplications, 1);
   varianceInLoopIndividual = zeros(numReplications, 1);
    
   
   
   
   % Create new vector to fill with MC replications
    
    parfor i=1:numReplications
        
        [y, x, u] = linearDynamicSimulateData(meanCoef+etaTrue/sqrt(T), sigmaSq, T, design, i); % draw data
        D = [y(T, 1); 1]; % conditional mean for unit with x_{T+1} = (-3, -3)
         betaScaled = meanCoef+etaTrue/sqrt(T);
         betaScaled =    dynamicChangeCoordinate1(lambdaRange(j), betaScaled);
        % Gradient block
 
        V = linearDynamicTrueCovariance(betaScaled(1,:), betaScaled(2,:),sigmaSq, N,design);    % compute population-fixed variances
       [psiMatrixTrue, ~, ~] = averagingPopulationPsi(etaTrue, V, D); % recompute population Psi
    %    psiMatrixTrue = averagingPopulationPsiOLD(etaTrue, V, D);
       [~, psiStein1] = averagingSteinWeights(etaTrue, compValue1, V, D); % create true matrix for Stein
       [~, psiStein2] = averagingSteinWeights(etaTrue, compValue2, V, D); % create true matrix for Stein
       [~, psiStein] = averagingSteinWeights(etaTrue, N, V, D); % create true matrix for Stein

        betaHat = linearStaticEstimators(y, x);% estimate coefficients
        etaEstimated = sqrt(T)*(betaHat - repmat(mean(betaHat,2),1, N));
        Vest = linearDynamicVarianceEstimator(y, x, betaHat); % estimate variance
        Dest= D; % plug into the gradient

        psiMatrixSample = averagingSamplePsi(betaHat*sqrt(T), Vest, Dest);% estimate sample psi
        [weightsStein1, psiSteinSample1] = averagingSteinWeights(betaHat*sqrt(T), compValue1, Vest, Dest);
        [weightsStein2, psiSteinSample2] = averagingSteinWeights(betaHat*sqrt(T), compValue2, Vest, Dest);
        [weightsStein, psiSteinSample] = averagingSteinWeights(betaHat*sqrt(T), N, Vest, Dest);
%       
%       
%         varianceInLoopOptimalP1(i) = averagingVariance(...
%                 averagingOptimalWeights(psiMatrixTrue(1:compValue1,1:compValue1), nW),...
%                 psiMatrixTrue(1:compValue1,1:compValue1)); 
%         varianceInLoopOptimalP2(i) = averagingVariance(...
%                 averagingOptimalWeights(psiMatrixTrue(1:compValue2,1:compValue2), nW),...
%                 psiMatrixTrue(1:compValue2,1:compValue2)); 
%         varianceInLoopOptimalPAll(i) = averagingVariance(...
%                 averagingOptimalWeights(psiMatrixTrue, nW), ...
%                 psiMatrixTrue);
            
        varianceInLoopSample1(i) = averagingVariance(...
                averagingOptimalWeights(psiMatrixSample(1:compValue1,1:compValue1), nW),...
                psiMatrixTrue(1:compValue1,1:compValue1)); 
        varianceInLoopSample2(i) = averagingVariance(...
                averagingOptimalWeights(psiMatrixSample(1:compValue2,1:compValue2), nW),...
                psiMatrixTrue(1:compValue2,1:compValue2)); 
        varianceInLoopSampleAll(i) = averagingVariance(...
                averagingOptimalWeights(psiMatrixSample, nW), psiMatrixTrue); 
        
          varianceInLoopEqual1(i) = averagingVariance(equalWeights1, psiMatrixTrue);
       varianceInLoopEqual2(i) = averagingVariance(equalWeights2, psiMatrixTrue);
       varianceInLoopEqualAll(i) = averagingVariance(equalWeightsAll, psiMatrixTrue);  

        varianceInLoopStein1(i) = averagingVariance(weightsStein1,...
                psiStein1); 
        varianceInLoopStein2(i) = averagingVariance(weightsStein2,...
                        psiStein2); 
        varianceInLoopStein(i) = averagingVariance(weightsStein,...
                        psiStein); 
        [aic1, mma1] =  linearAICWeights(betaHat,y, x, compValue1);
        [aic2, mma2] =  linearAICWeights(betaHat,y, x, compValue2);
        [aicAll, mmaAll] =  linearAICWeights(betaHat,y, x, N);
        varianceInLoopAIC1(i) = averagingVariance(aic1,...
                psiMatrixTrue(1:compValue1,1:compValue1)); 
        varianceInLoopAIC2(i) = averagingVariance(aic2,...
                psiMatrixTrue(1:compValue2,1:compValue2)); 
        varianceInLoopAIC(i) = averagingVariance(aicAll, psiMatrixTrue); 

        varianceInLoopMMA1(i) = averagingVariance(mma1,...
                psiMatrixTrue(1:compValue1,1:compValue1)); 
        varianceInLoopMMA2(i) = averagingVariance(mma2,...
                psiMatrixTrue(1:compValue2,1:compValue2)); 
        varianceInLoopMMA(i) = averagingVariance(mmaAll, psiMatrixTrue);
        varianceInLoopIndividual(i) = Dest'*V(:, :, 1)*Dest;    
            
        
        [j, i]
    end
     varianceSample1(j) = mean(varianceInLoopSample1);
     varianceSample2(j) = mean(varianceInLoopSample2);
     varianceSampleAll(j) = mean(varianceInLoopSampleAll);
     varianceStein1(j) =  mean(varianceInLoopStein1);
     varianceStein2(j) =  mean(varianceInLoopStein1);
     varianceStein(j) =  mean(varianceInLoopStein1);
     varianceAIC1(j) = mean(varianceInLoopAIC1);
     varianceAIC2(j) = mean(varianceInLoopAIC2);
     varianceAIC(j) = mean(varianceInLoopAIC);
     varianceMMA1(j) = mean(varianceInLoopMMA1);
     varianceMMA2(j) = mean(varianceInLoopMMA2);
     varianceMMA(j) = mean(varianceInLoopMMA);

     varianceIndividual(j) = mean(varianceInLoopIndividual);
         varianceEqual1(j) =mean(varianceInLoopEqual1);
     varianceEqual2(j) = mean(varianceInLoopEqual2);
     varianceEqualAll(j) = mean(varianceInLoopEqualAll);  
end

% Plot 1
figure


subplot(1,3,1)
    plot(lambdaRange,varianceSample1./varianceIndividual,'c','LineWidth',2) % optimal
    hold on
    plot(lambdaRange,varianceEqual1./varianceIndividual,'--')
    plot(lambdaRange, varianceStein1./varianceIndividual, '.-')
    plot(lambdaRange, varianceAIC1./varianceIndividual, '.-')
    plot(lambdaRange, varianceMMA1./varianceIndividual, '.-')
    legend('Plug-In', 'Equal','Stein-like', 'AIC','MMA', 'Location', 'northwest')
    title('N=5')
    xlim([min(lambdaRange), max(lambdaRange)])
    ylim([0, 2])
    xlabel('\lambda_1')

subplot(1,3,2)
    plot(lambdaRange,varianceSample2./varianceIndividual,'c','LineWidth',2) % optimal
    hold on
    plot(lambdaRange,varianceEqual2./varianceIndividual,'--')
    plot(lambdaRange, varianceStein2./varianceIndividual, '.-')
    plot(lambdaRange, varianceAIC2./varianceIndividual, '.-')
    plot(lambdaRange, varianceMMA2./varianceIndividual, '.-')
    legend('Plug-In', 'Equal','Stein-like', 'AIC','MMA', 'Location', 'northwest')
    title('N=25')
    xlim([min(lambdaRange), max(lambdaRange)])
    ylim([0, 2])
    xlabel('\lambda_1')

subplot(1,3,3)
    plot(lambdaRange,varianceSampleAll./varianceIndividual,'c','LineWidth',2) % optimal
    hold on
    plot(lambdaRange,varianceEqualAll./varianceIndividual,'--')
    plot(lambdaRange, varianceStein./varianceIndividual, '.-')
    plot(lambdaRange, varianceAIC./varianceIndividual, '.-')
    plot(lambdaRange, varianceMMA./varianceIndividual, '.-')
    legend('Plug-In', 'Equal','Stein-like', 'AIC','MMA', 'Location', 'northwest')
    title('N=100')
    xlim([min(lambdaRange), max(lambdaRange)])
    xlabel('\lambda_1')
    ylim([0, 2])

suptitle('Averaging Estiator, \mu(\theta_1) = E(y_{T+1}|y_T, x=1), ratio of MSE to individual estimator')
