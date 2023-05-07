

%% Define parameters
% Define parameters used in all simulations here

% Values for first and second plot
compValue1 = 5;
compValue2 = 25;

% Select a seed to drawing coefficients, 15
seedCoefficients =12133;

% For static model: select default mean and variance parameters for
% coefficients
meanBeta = [0;0];
varianceBeta = 1  ;
correlationBeta = -0.6;

% Simulation design choice
design =1;

% Data dimensions
N= 100;
T = 1000;

% Number of replications
numReplications = 100;

% Negative weights, set to 0 to enforce nonnegativity
nW = 0;

% Plot 1: MSE for changing the coefficient value

% Draw coefficients and variances
[betaOriginal,~ , sigmaSq] = linearStaticDrawCoefficients(N, T, meanBeta, varianceBeta,...
            correlationBeta, seedCoefficients);

betaOriginal  = -betaOriginal;
        
% Set range to be explored
beta1range = -5:0.05:5;
endOfRange = length(beta1range);

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


% Main loop: loop through parameter vector, draw multiple samples
for j=1:endOfRange
   
   % Change coordinate
   betaLoop =    staticChangeCoordinate1(beta1range(j), betaOriginal);

   % Gradient describes the parameter of interest
   D = [1; 0]; % first coefficient
%     D = [1; -2]; % conditional mean for unit with x_{T+1} = (1, -2)
%    D = [2*betaLoop(2,1); 2*betaLoop(1, 1)];  % norm of betad]

   etaTrue = betaLoop - repmat(meanBeta, 1, N); % obtain deviations from the mean
   V = linearStaticTrueCovariance(sigmaSq, N);    % compute population-fixed variances
%    [psiMatrixTrue, ~, ~] = averagingPopulationPsi(etaTrue, V, D); % recompute population Psi
   psiMatrixTrue = averagingPopulationPsiOLD(etaTrue, V, D);
   [~, psiStein1] = averagingSteinWeights(etaTrue, compValue1, V, D); % create true matrix for Stein
   [~, psiStein2] = averagingSteinWeights(etaTrue, compValue2, V, D); % create true matrix for Stein
   [~, psiStein] = averagingSteinWeights(etaTrue, N, V, D); % create true matrix for Stein

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
   varianceInLoopIndividual = zeros(numReplications, 1);
    
   
   
   
   % Create new vector to fill with MC replications
    
    parfor i=1:numReplications
        
        [y, x, u] = linearStaticSimulateData(meanBeta+ etaTrue/sqrt(T), sigmaSq, T, design, i); % draw data
        betaHat = linearStaticEstimators(y, x);% estimate coefficients
        Vest = linearStaticVarianceEstimator(y, x, betaHat); % estimate variance
        Dest= D; % plug into the gradient
%          Dest = [2*betaHat(1,1); 2*betaHat(2,1)];
        psiMatrixSample = averagingSamplePsi(betaHat*sqrt(T), Vest, Dest);% estimate sample psi
        [weightsStein1, psiSteinSample1] = averagingSteinWeights(betaHat*sqrt(T), compValue1, Vest, Dest);
        [weightsStein2, psiSteinSample2] = averagingSteinWeights(betaHat*sqrt(T), compValue2, Vest, Dest);
        [weightsStein, psiSteinSample] = averagingSteinWeights(betaHat*sqrt(T)  , N, Vest, Dest);
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
            
        
        
    end
    j
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
     varianceEqual1(j) = averagingVariance(equalWeights1, psiMatrixTrue);
     varianceEqual2(j) = averagingVariance(equalWeights2, psiMatrixTrue);
     varianceEqualAll(j) = averagingVariance(equalWeightsAll, psiMatrixTrue);    
end

% Plot 1
figure

subplot(1,3,1)
    plot(beta1range,varianceSample1./varianceIndividual,'c','LineWidth',2) % optimal
    hold on
    plot(beta1range,varianceEqual1./varianceIndividual,'--')
    plot(beta1range, varianceStein1./varianceIndividual, '.-')
    plot(beta1range, varianceAIC1./varianceIndividual, '.-')
    plot(beta1range, varianceMMA1./varianceIndividual, '.-')
    legend('Plug-In', 'Equal','Stein', 'AIC','MMA')
    title('N=5')
    xlim([min(beta1range), max(beta1range)])
    ylim([0, 2])
    xlabel('\eta')

subplot(1,3,2)
    plot(beta1range,varianceSample2./varianceIndividual,'c','LineWidth',2) % optimal
    hold on
    plot(beta1range,varianceEqual2./varianceIndividual,'--')
    plot(beta1range, varianceStein2./varianceIndividual, '.-')
    plot(beta1range, varianceAIC2./varianceIndividual, '.-')
    plot(beta1range, varianceMMA2./varianceIndividual, '.-')
    legend('Plug-In', 'Equal','Stein', 'AIC','MMA')
    title('N=25')
    xlim([min(beta1range), max(beta1range)])
    ylim([0, 2])
    xlabel('\eta')

subplot(1,3,3)
    plot(beta1range,varianceSampleAll./varianceIndividual,'c','LineWidth',2) % optimal
    hold on
    plot(beta1range,varianceEqualAll./varianceIndividual,'--')
    plot(beta1range, varianceStein./varianceIndividual, '.-')
    plot(beta1range, varianceAIC./varianceIndividual, '.-')
    plot(beta1range, varianceMMA./varianceIndividual, '.-')
    legend('Plug-In', 'Equal','Stein', 'AIC','MMA')
    title('N=100')
    xlim([min(beta1range), max(beta1range)])
    xlabel('\eta')
    ylim([0, 2])

% suptitle('Averaging Estimator, \mu(\beta) = \beta_1, ratio of MSE to individual estimator')
% suptitle('Averaging Estimator, \mu(\beta) = E(y|x_1=1, x_2=-2), ratio of MSE to individual estimator')
% suptitle('Averaging Estimator, \mu(\beta) = ||\beta||^2, ratio of MSE to individual estimator')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%% Plot 2: MSE for changing variance
% Create objects for later use
 % Create equal weights
 
        
sigma1range = 0.01:0.01:10;
endOfRange = length(sigma1range);
    

 equalWeights1 = zeros(N,1);
 equalWeights1(1:compValue1) = ones(compValue1,1)/5;
 equalWeights2 = zeros(N,1);
 equalWeights2(1:compValue2) = ones(compValue2,1)/compValue2;
 equalWeightsAll = ones(N,1)/(N);
 equalWeightsAll(1) =0;
 
 varianceIndividual = zeros(endOfRange, 1);
 varianceSample1 = zeros(endOfRange, 1);
 varianceSample2 = zeros(endOfRange, 1);
 varianceSampleAll = zeros(endOfRange, 1);
 varianceP1 = zeros(endOfRange, 1);
 varianceP2 = zeros(endOfRange, 1);
 variancePAll = zeros(endOfRange, 1);
 varianceStein1 = zeros(endOfRange, 1);
 varianceStein2 = zeros(endOfRange, 1);
 varianceStein = zeros(endOfRange, 1);
 varianceEqual1 = zeros(endOfRange, 1);
 varianceEqual2 = zeros(endOfRange, 1);
 varianceEqualAll = zeros(endOfRange, 1);
 varianceAIC1 = zeros(endOfRange, 1);
 varianceAIC2 = zeros(endOfRange, 1);
 varianceAIC = zeros(endOfRange, 1);
 varianceMMA1 = zeros(endOfRange, 1);
 varianceMMA2 = zeros(endOfRange, 1);
 varianceMMA = zeros(endOfRange, 1);

% Main loop

[betaOriginal, ~, sigmaSq] = linearStaticDrawCoefficients(N, T, meanBeta, varianceBeta,...
            correlationBeta, seedCoefficients);

for j=1:endOfRange
   

   betaLoop =  staticChangeScale(betaOriginal, sigma1range(j));
   % select new value for parameter
%    D = [1; 0];
   D = [1; -2]; % recompute the gradient D with new values
%    D = [2*betaLoop(2,1); 2*betaLoop(1, 1)];  % ratio

   etaTrue = betaLoop - repmat(meanBeta, 1, N); % obtain deviations from the mean
   V = linearStaticTrueCovariance(sigmaSq, N);    % compute population-fixed variances
   psiMatrixTrue = averagingPopulationPsiOLD(etaTrue , V, D); % recompute population Psi
   [~, psiStein1] = averagingSteinWeights(etaTrue, compValue1, V, D); % create true matrix for Stein
   [~, psiStein2] = averagingSteinWeights(etaTrue, compValue2, V, D); % create true matrix for Stein
   [~, psiStein] = averagingSteinWeights(etaTrue, N, V, D); % create true matrix for Stein

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
   varianceInLoopIndividual = zeros(numReplications, 1);
    
   
   
   
   % Create new vector to fill with MC replications
    
    parfor i=1:numReplications
        
        [y, x, u] = linearStaticSimulateData(betaLoop, sigmaSq, T, design, i); % draw data
        betaHat = linearStaticEstimators(y, x);% estimate coefficients
        Vest = linearStaticVarianceEstimator(y, x, betaHat); % estimate variance
        Dest= D; % plug into the gradient
%          Dest = [2*betaHat(1,1); 2*betaHat(2,1)];
        psiMatrixSample = averagingSamplePsi(betaHat, Vest, Dest);% estimate sample psi
        [weightsStein1, psiSteinSample1] = averagingSteinWeights(betaHat, compValue1, Vest, Dest);
        [weightsStein2, psiSteinSample2] = averagingSteinWeights(betaHat, compValue2, Vest, Dest);
        [weightsStein, psiSteinSample] = averagingSteinWeights(betaHat, N, Vest, Dest);
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
            
        
%         [j, i]
    end
    j
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
     varianceEqual1(j) = averagingVariance(equalWeights1, psiMatrixTrue);
     varianceEqual2(j) = averagingVariance(equalWeights2, psiMatrixTrue);
     varianceEqualAll(j) = averagingVariance(equalWeightsAll, psiMatrixTrue);    
end

% Plot 2
figure
subplot(1,3,1)

 plot(sigma1range,varianceSample1./varianceIndividual,'c','LineWidth',2) % optimal
 hold on
 plot(sigma1range,varianceEqual1./varianceIndividual,'--')
plot(sigma1range, varianceStein1./varianceIndividual, '.-')
 plot(sigma1range, varianceAIC1./varianceIndividual, '.-')
 plot(sigma1range, varianceMMA1./varianceIndividual, '.-')
legend('Plug-In', 'Equal','Stein', 'AIC','MMA')
title('N=5')
xlim([min(sigma1range), max(sigma1range)])
 ylim([0, 2])
xlabel('\sigma')

     subplot(1,3,2)
  plot(sigma1range,varianceSample2./varianceIndividual,'c','LineWidth',2) % optimal
hold on
  plot(sigma1range,varianceEqual2./varianceIndividual,'--')
plot(sigma1range, varianceStein2./varianceIndividual, '.-')
 plot(sigma1range, varianceAIC2./varianceIndividual, '.-')
 plot(sigma1range, varianceMMA2./varianceIndividual, '.-')
legend('Plug-In', 'Equal','Stein', 'AIC','MMA')
title('N=25')
xlim([min(sigma1range), max(sigma1range)])
 ylim([0, 2])
xlabel('\sigma')
subplot(1,3,3)

plot(sigma1range,varianceSampleAll./varianceIndividual,'c','LineWidth',2) % optimal
    hold on
    plot(sigma1range,varianceEqualAll./varianceIndividual,'--')
    plot(sigma1range, varianceStein./varianceIndividual, '.-')
    plot(sigma1range, varianceAIC./varianceIndividual, '.-')
    plot(sigma1range, varianceMMA./varianceIndividual, '.-')
    legend('Plug-In', 'Equal','Stein', 'AIC','MMA')
    title('N=100')
    xlim([min(sigma1range), max(sigma1range)])
    xlabel('\sigma')
    ylim([0, 2])

% suptitle('Averaging Estimator, \mu(\beta) = \beta_1, ratio of MSE to individual estimator')
suptitle('Averaging Estimator, \mu(\beta) = E(y|x_1=1, x_2=-2), ratio of MSE to individual estimator')
% suptitle('Averaging Estimator, \mu(\beta) = ||\beta||^2, ratio of MSE to individual estimator')
