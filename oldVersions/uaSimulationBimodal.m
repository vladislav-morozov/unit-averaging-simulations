
%%

subplot(1,3,1)
    plot(lambda1range,mseSample1./mseIndividual,'c','LineWidth',2) % optimal
    hold on
    plot(lambda1range,mseEqual1./mseIndividual,'--')
    plot(lambda1range, mseStein1./mseIndividual, '.-')
    plot(lambda1range, mseAIC1./mseIndividual, '.-')
    plot(lambda1range, mseMMA1./mseIndividual, '.-')
    legend('Plug-In', 'Equal','Stein-like', 'AIC','MMA', 'Location', 'northeast')
    title('N=5')
    xlim([min(lambda1range), max(lambda1range)])
    ylim([0, 2])
    xlabel('\eta_1')

subplot(1,3,2)
    plot(lambda1range,mseSample2./mseIndividual,'c','LineWidth',2) % optimal
    hold on
    plot(lambda1range,mseEqual2./mseIndividual,'--')
    plot(lambda1range, mseStein2./mseIndividual, '.-')
    plot(lambda1range, mseAIC2./mseIndividual, '.-')
    plot(lambda1range, mseMMA2./mseIndividual, '.-')
    legend('Plug-In', 'Equal','Stein-like', 'AIC','MMA', 'Location', 'northeast')
    title('N=25')
    xlim([min(lambda1range), max(lambda1range)])
    ylim([0, 2])
    xlabel('\eta_1')

subplot(1,3,3)
    plot(lambda1range,mseSampleAll./mseIndividual,'c','LineWidth',2) % optimal
    hold on
    plot(lambda1range,mseEqualAll./mseIndividual,'--')
    plot(lambda1range, mseStein./mseIndividual, '.-')
    plot(lambda1range, mseAIC./mseIndividual, '.-')
    plot(lambda1range, mseMMA./mseIndividual, '.-')
    legend('Plug-In', 'Equal','Stein-like', 'AIC','MMA', 'Location', 'northeast')
    title('N=100')
    xlim([min(lambda1range), max(lambda1range)])
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
[thetaSample,~ , sigmaSq] = linearDynamicDrawCoefficients(N, T, design,... 
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
 mseIndividual = zeros(endOfRange, 1);
    % Plug-in optimized
 mseSample1 = zeros(endOfRange, 1);
 mseSample2 = zeros(endOfRange, 1);
 mseSampleAll = zeros(endOfRange, 1);
    % Population optimal
 mseP1 = zeros(endOfRange, 1);
 mseP2 = zeros(endOfRange, 1);
 msePAll = zeros(endOfRange, 1);
    % Stein-type
 mseStein1 = zeros(endOfRange, 1);
 mseStein2 = zeros(endOfRange, 1);
 mseStein = zeros(endOfRange, 1);
    % Equal
 mseEqual1 = zeros(endOfRange, 1);
 mseEqual2 = zeros(endOfRange, 1);
 mseEqualAll = zeros(endOfRange, 1);
    % AIC/BIC
 mseAIC1 = zeros(endOfRange, 1);
 mseAIC2 = zeros(endOfRange, 1);
 mseAIC = zeros(endOfRange, 1);
    % MMA
 mseMMA1 = zeros(endOfRange, 1);
 mseMMA2 = zeros(endOfRange, 1);
 mseMMA = zeros(endOfRange, 1);


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
       betaLoop =    thetaSample;
       betaLoop(2,:) = linearDynamicDrawBimodal(a,-a, 1, 1, 0.5, N,i);
       etaTrue = betaLoop - repmat(meanCoef, 1, N); % obtain deviations from the mean

       % Rescale coefficients and set the value of the first to the value explored
       thetaScaled = meanCoef+etaTrue/sqrt(T);
    %    betaScaled =    dynamicChangeCoordinate1(lambdaRange(j), betaScaled);

       % Choose gradient, gradient describes which parameter is estimated
%        D = [0; 1] % beta_1
%        D = [betaScaled(2, 1)/(1-betaScaled(1, 1))^2; 1/(1-betaScaled(1, 1))]; % long-run
%        D = [2*betaScaled(2,1); 2*betaScaled(1, 1)] % norm
           

   
        % Draw data, estimate coefficients and the gradent
        [y, x, u] = linearDynamicSimulateData(meanCoef+etaTrue/sqrt(T), sigmaSq, T, design, i); % draw data
       
        d0 = [y(T, 1); 1]; % forecast
        
        % Compute variance and population Psi matrices
       V = linearDynamicTrueCovariance(thetaScaled(1,:), thetaScaled(2,:),sigmaSq, N,design);    % compute population-fixed variances
       [psiMatrixTrue, ~, ~] = averagingPopulationPsi(etaTrue, V, d0); % recompute population Psi
       % Special computations for Stein-type estimator 
       [~, psiStein1] = averagingSteinWeights(etaTrue, compValue1, V, d0); % create true matrix for Stein
       [~, psiStein2] = averagingSteinWeights(etaTrue, compValue2, V, d0); % create true matrix for Stein
       [~, psiStein] = averagingSteinWeights(etaTrue, N, V, d0); % create true matrix for Stein

        
        
        thetaHat = linearStaticEstimators(y, x);% estimate coefficients
        etaEst = sqrt(T)*(thetaHat - repmat(mean(thetaHat,2),1, N));
        Vest = linearDynamicVarianceEstimator(y, x, thetaHat); % estimate variance
        Dest= d0; % plug into the gradient
%         Dest = [2*betaHat(1,1); 2*betaHat(2,1)];
%         Dest = [betaHat(2, 1)/(1-betaHat(1, 1))^2; 1/(1-betaHat(1, 1))];

        % Compute Psi matrics based on data
        psiMatrixSample = averagingSamplePsi(thetaHat*sqrt(T), Vest, Dest);% estimate sample psi
        [weightsStein1, psiSteinSample1] = averagingSteinWeights(thetaHat*sqrt(T), compValue1, Vest, Dest);
        [weightsStein2, psiSteinSample2] = averagingSteinWeights(thetaHat*sqrt(T), compValue2, Vest, Dest);
        [weightsStein, psiSteinSample] = averagingSteinWeights(thetaHat*sqrt(T), N, Vest, Dest);

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
        [aic1, mma1] =  linearAICWeights(thetaHat,y, x, compValue1);
        [aic2, mma2] =  linearAICWeights(thetaHat,y, x, compValue2);
        [aicAll, mmaAll] =  linearAICWeights(thetaHat,y, x, N);
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
     mseSample1(j) = mean(varianceInLoopSample1);
     mseSample2(j) = mean(varianceInLoopSample2);
     mseSampleAll(j) = mean(varianceInLoopSampleAll);
     mseStein1(j) =  mean(varianceInLoopStein1);
     mseStein2(j) =  mean(varianceInLoopStein1);
     mseStein(j) =  mean(varianceInLoopStein1);
     mseAIC1(j) = mean(varianceInLoopAIC1);
     mseAIC2(j) = mean(varianceInLoopAIC2);
     mseAIC(j) = mean(varianceInLoopAIC);
     mseMMA1(j) = mean(varianceInLoopMMA1);
     mseMMA2(j) = mean(varianceInLoopMMA2);
     mseMMA(j) = mean(varianceInLoopMMA);
     mseP1(j) = averagingVariance(...
                    averagingOptimalWeights(psiMatrixTrue(1:compValue1,1:compValue1), nW),...
                 psiMatrixTrue(1:compValue1,1:compValue1));
     mseP2(j) = averagingVariance(...
                    averagingOptimalWeights(psiMatrixTrue(1:compValue2,1:compValue2), nW),...
                  psiMatrixTrue(1:compValue2,1:compValue2));
     msePAll(j) = averagingVariance(...
                    averagingOptimalWeights(psiMatrixTrue, nW),...
                 psiMatrixTrue);
     mseIndividual(j) = mean(varianceInLoopIndividual);

          mseIndividual(j) = mean(varianceInLoopIndividual);
         mseEqual1(j) =mean(varianceInLoopEqual1);
     mseEqual2(j) = mean(varianceInLoopEqual2);
     mseEqualAll(j) = mean(varianceInLoopEqualAll);  
end

% Plot 3
figure

subplot(1,3,1)
    plot(aRange,mseSample1./mseIndividual,'c','LineWidth',2) % optimal
    hold on
    plot(aRange,mseEqual1./mseIndividual,'--')
    plot(aRange, mseStein1./mseIndividual, '.-')
    plot(aRange, mseAIC1./mseIndividual, '.-')
    plot(aRange, mseMMA1./mseIndividual, '.-')
    legend('Plug-In', 'Equal','Stein-like', 'AIC','MMA', 'Location', 'northeast')
    title('N=5')
    xlim([min(aRange), max(aRange)])
    ylim([0, 2])
    xlabel('a')

subplot(1,3,2)
    plot(aRange,mseSample2./mseIndividual,'c','LineWidth',2) % optimal
    hold on
    plot(aRange,mseEqual2./mseIndividual,'--')
    plot(aRange, mseStein2./mseIndividual, '.-')
    plot(aRange, mseAIC2./mseIndividual, '.-')
    plot(aRange, mseMMA2./mseIndividual, '.-')
    legend('Plug-In', 'Equal','Stein-like', 'AIC','MMA', 'Location', 'northeast')
    title('N=25')
    xlim([min(aRange), max(aRange)])
    ylim([0, 2])
    xlabel('a')

subplot(1,3,3)
    plot(aRange,mseSampleAll./mseIndividual,'c','LineWidth',2) % optimal
    hold on
    plot(aRange,mseEqualAll./mseIndividual,'--')
    plot(aRange, mseStein./mseIndividual, '.-')
    plot(aRange, mseAIC./mseIndividual, '.-')
    plot(aRange, mseMMA./mseIndividual, '.-')
    legend('Plug-In', 'Equal','Stein-like', 'AIC','MMA', 'Location', 'northeast')
    title('N=100')
    xlim([min(aRange), max(aRange)])
    xlabel('a')
    ylim([0, 2])
% 
% suptitle('Averaging Estimator, \mu(\theta_1) = \beta_1, ratio of MSE to individual estimator')
suptitle('Averaging Estiator, \mu(\theta_1) = E(y_{T+1}|y_T=1=-1, x=1), ratio of MSE to individual estimator')
%  suptitle('Averaging Estimator, \mu(\theta_1) = {\beta_1}/{1-\lambda_1}, ratio of MSE to individual estimator')
