%% Plot 1: MSE for changing the coefficient value

% Draw coefficients and variances
[betaOriginal,~ , sigmaSq] = linearDynamicDrawCoefficients(N, T, design,... 
            meanBeta, varianceBeta,   seedCoefficients);

% Set range to be explored
lambdaRange = 0:0.01:0.99;
endOfRange = length(lambdaRange);

% Create objects for later use
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


% Main loop: loop through parameter vector, draw multiple samples for each
% value

for j=1:endOfRange
   
   % Change coordinate and obtain deviations
   betaLoop =    dynamicChangeCoordinate1(lambdaRange(j), betaOriginal);
   etaTrue = betaLoop - repmat(meanCoef, 1, N); % obtain deviations from the mean
   
   % Rescale coefficients and set the value of the first to the value explored
   betaScaled = meanCoef+etaTrue/sqrt(T);
%    betaScaled =    dynamicChangeCoordinate1(lambdaRange(j), betaScaled);
   
   % Choose gradient, gradient describes which parameter is estimated
   % D = [1; 0] % lambda_1
   D = [betaScaled(2, 1)/(1-betaScaled(1, 1))^2; 1/(1-betaScaled(1, 1))]; % long-run
   % D = [2*betaScaled(2,1); 2*betaScaled(1, 1)] % norm

   % Compute variance and population Psi matrices
   V = linearDynamicTrueCovariance(betaScaled(1,:), betaScaled(2,:),sigmaSq, N,design);    % compute population-fixed variances
   [psiMatrixTrue, ~, ~] = averagingPopulationPsi(etaTrue, V, D); % recompute population Psi
   % Special computations for Stein-type estimator 
   [~, psiStein1] = averagingSteinWeights(etaTrue, compValue1, V, D); % create true matrix for Stein
   [~, psiStein2] = averagingSteinWeights(etaTrue, compValue2, V, D); % create true matrix for Stein
   [~, psiStein] = averagingSteinWeights(etaTrue, N, V, D); % create true matrix for Stein

   
   % Recreate temporary variance vectors
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

   % Inner loop: drawing data
    parfor i=1:numReplications
        
        % Draw data, estimate coefficients and the gradent
        [y, x, u] = linearDynamicSimulateData(meanCoef+etaTrue/sqrt(T), sigmaSq, T, design, i); % draw data
        betaHat = linearStaticEstimators(y, x);% estimate coefficients
        etaEstimated = sqrt(T)*(betaHat - repmat(mean(betaHat,2),1, N));
        Vest = linearDynamicVarianceEstimator(y, x, betaHat); % estimate variance
%         Dest= D; % plug into the gradient
%         Dest = [2*betaHat(1,1); 2*betaHat(2,1)];
        Dest = [betaHat(2, 1)/(1-betaHat(1, 1))^2; 1/(1-betaHat(1, 1))];

        % Compute Psi matrics based on data
        psiMatrixSample = averagingSamplePsi(betaHat*sqrt(T), Vest, Dest);% estimate sample psi
        [weightsStein1, psiSteinSample1] = averagingSteinWeights(betaHat*sqrt(T), compValue1, Vest, Dest);
        [weightsStein2, psiSteinSample2] = averagingSteinWeights(betaHat*sqrt(T), compValue2, Vest, Dest);
        [weightsStein, psiSteinSample] = averagingSteinWeights(betaHat*sqrt(T), N, Vest, Dest);

        % Compute all variances
        varianceInLoopSample1(i) = averagingVariance(...
                averagingOptimalWeights(psiMatrixSample(1:compValue1,1:compValue1), nW),...
                psiMatrixTrue(1:compValue1,1:compValue1)); 
        varianceInLoopSample2(i) = averagingVariance(...
                averagingOptimalWeights(psiMatrixSample(1:compValue2,1:compValue2), nW),...
                psiMatrixTrue(1:compValue2,1:compValue2)); 
        varianceInLoopSampleAll(i) = averagingVariance(...
                averagingOptimalWeights(psiMatrixSample, nW), psiMatrixTrue); 
        

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
        
       varianceInLoopEqual1(i) = averagingVariance(equalWeights1, psiMatrixTrue);
       varianceInLoopEqual2(i) = averagingVariance(equalWeights2, psiMatrixTrue);
       varianceInLoopEqualAll(i) = averagingVariance(equalWeightsAll, psiMatrixTrue);  

        
    end
    % Average over data sets for each coordinate
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
     varianceP1(j) = averagingVariance(...
                    averagingOptimalWeights(psiMatrixTrue(1:compValue1,1:compValue1), nW),...
                 psiMatrixTrue(1:compValue1,1:compValue1));
     varianceP2(j) = averagingVariance(...
                    averagingOptimalWeights(psiMatrixTrue(1:compValue2,1:compValue2), nW),...
                  psiMatrixTrue(1:compValue2,1:compValue2));
     variancePAll(j) = averagingVariance(...
                    averagingOptimalWeights(psiMatrixTrue, nW),...
                 psiMatrixTrue);
     varianceIndividual(j) = mean(varianceInLoopIndividual);

          varianceIndividual(j) = mean(varianceInLoopIndividual);
         varianceEqual1(j) =mean(varianceInLoopEqual1);
     varianceEqual2(j) = mean(varianceInLoopEqual2);
     varianceEqualAll(j) = mean(varianceInLoopEqualAll);  
end

%% Plot 1
figure

subplot(1,3,1)
    plot(lambdaRange,varianceSample1./varianceIndividual,'c','LineWidth',2) % optimal
    hold on
    plot(lambdaRange,varianceEqual1./varianceIndividual,'--')
    plot(lambdaRange, varianceStein1./varianceIndividual, '.-')
    plot(lambdaRange, varianceAIC1./varianceIndividual, '.-')
    plot(lambdaRange, varianceMMA1./varianceIndividual, '.-')
    legend('Plug-In', 'Equal','Stein-like', 'AIC','MMA', 'Location', 'northeast')
    title('N=5')
    xlim([min(lambdaRange), max(lambdaRange)])
    ylim([0, 2])
    xlabel('\eta_1')

subplot(1,3,2)
    plot(lambdaRange,varianceSample2./varianceIndividual,'c','LineWidth',2) % optimal
    hold on
    plot(lambdaRange,varianceEqual2./varianceIndividual,'--')
    plot(lambdaRange, varianceStein2./varianceIndividual, '.-')
    plot(lambdaRange, varianceAIC2./varianceIndividual, '.-')
    plot(lambdaRange, varianceMMA2./varianceIndividual, '.-')
    legend('Plug-In', 'Equal','Stein-like', 'AIC','MMA', 'Location', 'northeast')
    title('N=25')
    xlim([min(lambdaRange), max(lambdaRange)])
    ylim([0, 2])
    xlabel('\eta_1')

subplot(1,3,3)
    plot(lambdaRange,varianceSampleAll./varianceIndividual,'c','LineWidth',2) % optimal
    hold on
    plot(lambdaRange,varianceEqualAll./varianceIndividual,'--')
    plot(lambdaRange, varianceStein./varianceIndividual, '.-')
    plot(lambdaRange, varianceAIC./varianceIndividual, '.-')
    plot(lambdaRange, varianceMMA./varianceIndividual, '.-')
    legend('Plug-In', 'Equal','Stein-like', 'AIC','MMA', 'Location', 'northeast')
    title('N=100')
    xlim([min(lambdaRange), max(lambdaRange)])
    xlabel('\eta_1')
    ylim([0, 2])
% 
suptitle('Averaging Estimator, \mu(\theta_1) = \lambda_1, ratio of MSE to individual estimator')
% suptitle('Averaging Estiator, \mu(\theta_1) = E(y_{T+1}|y_T=1=-1, x=1), ratio of MSE to individual estimator')
%  suptitle('Averaging Estimator, \mu(\theta_1) = {\beta_1}/{1-\lambda_1}, ratio of MSE to individual estimator')


 %% For one-step ahead forecasts, run this instead
 averagingLinearDynamicMConeStepAheadlambda
 
%% Effect of bimodality on estimation


% Draw coefficients and variances
[betaOriginal,~ , sigmaSq] = linearDynamicDrawCoefficients(N, T, design,... 
            meanBeta, varianceBeta,   seedCoefficients);

% Set range to be explored
aRange = 0:0.125:10;
endOfRange = length(aRange);

% Create objects for later use
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


% Main loop: loop through parameter vector, draw multiple samples for each
% value

for j=1:endOfRange
   

   a = aRange(j);
   % Recreate temporary variance vectors
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

   % Inner loop: drawing data
    parfor i=1:numReplications
        
           % Change coordinate and obtain deviations
       betaLoop =    betaOriginal;
       betaLoop(2,:) = linearDynamicDrawBimodal(a,-a, 1, 1, 0.5, N,i);
       etaTrue = betaLoop - repmat(meanCoef, 1, N); % obtain deviations from the mean

       % Rescale coefficients and set the value of the first to the value explored
       betaScaled = meanCoef+etaTrue/sqrt(T);
    %    betaScaled =    dynamicChangeCoordinate1(lambdaRange(j), betaScaled);

       % Choose gradient, gradient describes which parameter is estimated
%        D = [0; 1] % beta_1
%        D = [betaScaled(2, 1)/(1-betaScaled(1, 1))^2; 1/(1-betaScaled(1, 1))]; % long-run
%        D = [2*betaScaled(2,1); 2*betaScaled(1, 1)] % norm
           

   
        % Draw data, estimate coefficients and the gradent
        [y, x, u] = linearDynamicSimulateData(meanCoef+etaTrue/sqrt(T), sigmaSq, T, design, i); % draw data
       
        D = [y(T, 1); 1]; % forecast
        
        % Compute variance and population Psi matrices
       V = linearDynamicTrueCovariance(betaScaled(1,:), betaScaled(2,:),sigmaSq, N,design);    % compute population-fixed variances
       [psiMatrixTrue, ~, ~] = averagingPopulationPsi(etaTrue, V, D); % recompute population Psi
       % Special computations for Stein-type estimator 
       [~, psiStein1] = averagingSteinWeights(etaTrue, compValue1, V, D); % create true matrix for Stein
       [~, psiStein2] = averagingSteinWeights(etaTrue, compValue2, V, D); % create true matrix for Stein
       [~, psiStein] = averagingSteinWeights(etaTrue, N, V, D); % create true matrix for Stein

        
        
        betaHat = linearStaticEstimators(y, x);% estimate coefficients
        etaEstimated = sqrt(T)*(betaHat - repmat(mean(betaHat,2),1, N));
        Vest = linearDynamicVarianceEstimator(y, x, betaHat); % estimate variance
        Dest= D; % plug into the gradient
%         Dest = [2*betaHat(1,1); 2*betaHat(2,1)];
%         Dest = [betaHat(2, 1)/(1-betaHat(1, 1))^2; 1/(1-betaHat(1, 1))];

        % Compute Psi matrics based on data
        psiMatrixSample = averagingSamplePsi(betaHat*sqrt(T), Vest, Dest);% estimate sample psi
        [weightsStein1, psiSteinSample1] = averagingSteinWeights(betaHat*sqrt(T), compValue1, Vest, Dest);
        [weightsStein2, psiSteinSample2] = averagingSteinWeights(betaHat*sqrt(T), compValue2, Vest, Dest);
        [weightsStein, psiSteinSample] = averagingSteinWeights(betaHat*sqrt(T), N, Vest, Dest);

        % Compute all variances
        varianceInLoopSample1(i) = averagingVariance(...
                averagingOptimalWeights(psiMatrixSample(1:compValue1,1:compValue1), nW),...
                psiMatrixTrue(1:compValue1,1:compValue1)); 
        varianceInLoopSample2(i) = averagingVariance(...
                averagingOptimalWeights(psiMatrixSample(1:compValue2,1:compValue2), nW),...
                psiMatrixTrue(1:compValue2,1:compValue2)); 
        varianceInLoopSampleAll(i) = averagingVariance(...
                averagingOptimalWeights(psiMatrixSample, nW), psiMatrixTrue); 
        

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
        
       varianceInLoopEqual1(i) = averagingVariance(equalWeights1, psiMatrixTrue);
       varianceInLoopEqual2(i) = averagingVariance(equalWeights2, psiMatrixTrue);
       varianceInLoopEqualAll(i) = averagingVariance(equalWeightsAll, psiMatrixTrue);  

       [j, i] 
    end
    % Average over data sets for each coordinate
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
     varianceP1(j) = averagingVariance(...
                    averagingOptimalWeights(psiMatrixTrue(1:compValue1,1:compValue1), nW),...
                 psiMatrixTrue(1:compValue1,1:compValue1));
     varianceP2(j) = averagingVariance(...
                    averagingOptimalWeights(psiMatrixTrue(1:compValue2,1:compValue2), nW),...
                  psiMatrixTrue(1:compValue2,1:compValue2));
     variancePAll(j) = averagingVariance(...
                    averagingOptimalWeights(psiMatrixTrue, nW),...
                 psiMatrixTrue);
     varianceIndividual(j) = mean(varianceInLoopIndividual);

          varianceIndividual(j) = mean(varianceInLoopIndividual);
         varianceEqual1(j) =mean(varianceInLoopEqual1);
     varianceEqual2(j) = mean(varianceInLoopEqual2);
     varianceEqualAll(j) = mean(varianceInLoopEqualAll);  
end

% Plot 3
figure

subplot(1,3,1)
    plot(aRange,varianceSample1./varianceIndividual,'c','LineWidth',2) % optimal
    hold on
    plot(aRange,varianceEqual1./varianceIndividual,'--')
    plot(aRange, varianceStein1./varianceIndividual, '.-')
    plot(aRange, varianceAIC1./varianceIndividual, '.-')
    plot(aRange, varianceMMA1./varianceIndividual, '.-')
    legend('Plug-In', 'Equal','Stein-like', 'AIC','MMA', 'Location', 'northeast')
    title('N=5')
    xlim([min(aRange), max(aRange)])
    ylim([0, 2])
    xlabel('a')

subplot(1,3,2)
    plot(aRange,varianceSample2./varianceIndividual,'c','LineWidth',2) % optimal
    hold on
    plot(aRange,varianceEqual2./varianceIndividual,'--')
    plot(aRange, varianceStein2./varianceIndividual, '.-')
    plot(aRange, varianceAIC2./varianceIndividual, '.-')
    plot(aRange, varianceMMA2./varianceIndividual, '.-')
    legend('Plug-In', 'Equal','Stein-like', 'AIC','MMA', 'Location', 'northeast')
    title('N=25')
    xlim([min(aRange), max(aRange)])
    ylim([0, 2])
    xlabel('a')

subplot(1,3,3)
    plot(aRange,varianceSampleAll./varianceIndividual,'c','LineWidth',2) % optimal
    hold on
    plot(aRange,varianceEqualAll./varianceIndividual,'--')
    plot(aRange, varianceStein./varianceIndividual, '.-')
    plot(aRange, varianceAIC./varianceIndividual, '.-')
    plot(aRange, varianceMMA./varianceIndividual, '.-')
    legend('Plug-In', 'Equal','Stein-like', 'AIC','MMA', 'Location', 'northeast')
    title('N=100')
    xlim([min(aRange), max(aRange)])
    xlabel('a')
    ylim([0, 2])
% 
% suptitle('Averaging Estimator, \mu(\theta_1) = \beta_1, ratio of MSE to individual estimator')
suptitle('Averaging Estiator, \mu(\theta_1) = E(y_{T+1}|y_T=1=-1, x=1), ratio of MSE to individual estimator')
%  suptitle('Averaging Estimator, \mu(\theta_1) = {\beta_1}/{1-\lambda_1}, ratio of MSE to individual estimator')
