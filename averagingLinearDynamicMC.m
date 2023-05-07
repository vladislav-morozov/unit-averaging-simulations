

%% Define parameters
% Define parameters used in all simulations here

% Values for first and second plot
compValue1 = 5;
compValue2 = 25;

% Select a seed to drawing coefficients, 15
seedCoefficients = 15;
% For Dynamic model: select default mean and variance parameters for
% coefficients
meanBeta = 1;
varianceBeta = 1;
correlationBeta = 0.6;
meanCoef = [0; meanBeta];

% Simulation design choice
design =2;

% Data dimensions
N= 100;
T = 100;

% Number of replications
numReplications = 200;

% Negative weights, set to 0 to enforce nonnegativity
nW = 0;

% Plot 1: MSE for changing the coefficient value

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


% Main loop: loop through parameter vector, draw multiple samples
for j=1:endOfRange
   
   % Change coordinate
   betaLoop =    dynamicChangeCoordinate1(lambdaRange(j), betaOriginal);
    
   % Gradient describes the parameter of interest
%    D = [1; 0]; % first coefficient
    D=[0; 1];
%     D = [betaLoop(2, 1)/(1-betaLoop(1, 1))^2; 1/(1-betaLoop(1, 1))]; % LR
%     D = [2*betaLoop(2,1); 2*betaLoop(1, 1)];  % norm of beta
 
   etaTrue = betaLoop - repmat(meanCoef, 1, N); % obtain deviations from the mean
   V = linearDynamicTrueCovariance(betaLoop(1,:), betaLoop(2,:),sigmaSq, N,design);    % compute population-fixed variances
   [psiMatrixTrue, ~, ~] = averagingPopulationPsi(etaTrue, V, D); % recompute population Psi
%    psiMatrixTrue = averagingPopulationPsiOLD(etaTrue, V, D);
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
        
        [y, x, u] = linearDynamicSimulateData(meanCoef+etaTrue/sqrt(T), sigmaSq, T, design, i); % draw data
        betaHat = linearStaticEstimators(y, x);% estimate coefficients
        etaEstimated = sqrt(T)*(betaHat - repmat(mean(betaHat,2),1, N));
        Vest = linearDynamicVarianceEstimator(y, x, betaHat); % estimate variance
%         Dest= D; % plug into the gradient
%          Dest = [2*betaHat(1,1); 2*betaHat(2,1)];
        Dest = [betaHat(2, 1)/(1-betaHat(1, 1))^2; 1/(1-betaHat(1, 1))];
%     
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
    xlabel('\lambda_1')

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
    xlabel('\lambda_1')

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
    xlabel('\lambda_1')
    ylim([0, 2])
% 
% suptitle('Averaging Estimator, \mu(\theta_1) = \lambda_1, ratio of MSE to individual estimator')
% suptitle('Averaging Estiator, \mu(\theta_1) = E(y_{T+1}|y_T=1=-1, x=1), ratio of MSE to individual estimator')
 suptitle('Averaging Estimator, \mu(\theta_1) = {\beta_1}/{1-\lambda_1}, ratio of MSE to individual estimator')

%% Call a subfile with a slightly slower loop for forecast
averagingLinearDynamicMConeStepAhead

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%% Plot 2: Averaging over distribution of eta as well
% Create objects for later use
 % Create equal weights
 
meanBeta = 1;
varianceBeta = 1;
correlationBeta = 0.8;
meanCoef = [0; meanBeta];

% Simulation design choice
design =1;

% Replication numbers
numDrawEta = 100;
numDrawData = 200;

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

[betaOriginal, ~, sigmaSq] = linearDynamicDrawCoefficients(N, T, meanBeta, varianceBeta,...
            correlationBeta, seedCoefficients);

% Outer loop: go through parameter space explored        
for j=1:endOfRange

     
     varianceMIndividual = zeros(numDrawEta, 1);
     varianceMSample1 = zeros(numDrawEta, 1);
     varianceMSample2 = zeros(numDrawEta, 1);
     varianceMSampleAll = zeros(numDrawEta, 1);
     varianceMP1 = zeros(numDrawEta, 1);
     varianceMP2 = zeros(numDrawEta, 1);
     varianceMPAll = zeros(numDrawEta, 1);
     varianceMStein1 = zeros(numDrawEta, 1);
     varianceMStein2 = zeros(numDrawEta, 1);
     varianceMStein = zeros(numDrawEta, 1);
     varianceMEqual1 = zeros(numDrawEta, 1);
     varianceMEqual2 = zeros(numDrawEta, 1);
     varianceMEqualAll = zeros(numDrawEta, 1);
     varianceMAIC1 = zeros(numDrawEta, 1);
     varianceMAIC2 = zeros(numDrawEta, 1);
     varianceMAIC = zeros(numDrawEta, 1);
     varianceMMMA1 = zeros(numDrawEta, 1);
     varianceMMMA2 = zeros(numDrawEta, 1);
     varianceMMMA = zeros(numDrawEta, 1);
         % Middle loop: draw coefficients and produce average variances for
    % given coefficient set
    for k=1:numDrawEta
       [betaOriginal, ~, sigmaSq] = linearDynamicDrawCoefficients(N, T, meanBeta, varianceBeta,...
                correlationBeta, k);
       betaLoop =    dynamicChangeCoordinate1(lambdaRange(j), betaOriginal);
       % select new value for parameter
       D = [1; 0];
        
       etaTrue = betaLoop - repmat(meanCoef, 1, N); % obtain deviations from the mean
       V = linearDynamicTrueCovariance(betaLoop(1,:), betaLoop(2,:),sigmaSq, N, design);    % compute population-fixed variances
       [psiMatrixTrue, ~, ~] = averagingPopulationPsi(etaTrue, V, D); % recompute population Psi
    %    psiMatrixTrue = averagingPopulationPsiOLD(etaTrue, V, D);
       [~, psiStein1] = averagingSteinWeights(etaTrue, compValue1, V, D); % create true matrix for Stein
       [~, psiStein2] = averagingSteinWeights(etaTrue, compValue2, V, D); % create true matrix for Stein
       [~, psiStein] = averagingSteinWeights(etaTrue, N, V, D); % create true matrix for Stein

       % Recreate variance vectors
       varianceInLoopSample1 = zeros(numDrawData, 1);
       varianceInLoopSample2 = zeros(numDrawData, 1);
       varianceInLoopSampleAll = zeros(numDrawData, 1);
       varianceInLoopStein1 = zeros(numDrawData, 1);
       varianceInLoopAIC1 = zeros(numDrawData, 1);
       varianceInLoopAIC2 = zeros(numDrawData, 1);
       varianceInLoopAIC = zeros(numDrawData, 1);
       varianceInLoopMMA1 = zeros(numDrawData, 1);
       varianceInLoopMMA2 = zeros(numDrawData, 1);
       varianceInLoopMMA = zeros(numDrawData, 1);
       varianceInLoopOptimalP1 = zeros(numDrawData, 1);
       varianceInLoopOptimalP2 = zeros(numDrawData, 1);
       varianceInLoopOptimalPAll = zeros(numDrawData, 1);
       varianceInLoopIndividual = zeros(numDrawData, 1);

       
        % Internal loop: draw datasets for coefficients
        parfor i=1:numDrawData

            [y, x, u] = linearDynamicSimulateData(meanCoef+etaTrue/sqrt(T), sigmaSq, T, design, i); % draw data
            betaHat = linearStaticEstimators(y, x);% estimate coefficients
            etaEstimated = sqrt(T)*(betaHat - repmat(mean(betaHat,2),1, N));
            Vest = linearDynamicVarianceEstimator(y, x, betaHat); % estimate variance
            Dest= D; % plug into the gradient
    %          Dest = [2*betaHat(1,1); 2*betaHat(2,1)];
%             Dest = [betaHat(2, 1)/(1-betaHat(1, 1))^2; 1/(1-betaHat(1, 1))];
    %     
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



            [j,k, i]
        end
    
         varianceMSample1(k) = mean(varianceInLoopSample1);
         varianceMSample2(k) = mean(varianceInLoopSample2);
         varianceMSampleAll(k) = mean(varianceInLoopSampleAll);
         varianceMStein1(k) =  mean(varianceInLoopStein1);
         varianceMStein2(k) =  mean(varianceInLoopStein1);
         varianceMStein(k) =  mean(varianceInLoopStein1);
         varianceMAIC1(k) = mean(varianceInLoopAIC1);
         varianceMAIC2(k) = mean(varianceInLoopAIC2);
         varianceMAIC(k) = mean(varianceInLoopAIC);
         varianceMMMA1(k) = mean(varianceInLoopMMA1);
         varianceMMMA2(k) = mean(varianceInLoopMMA2);
         varianceMMMA(k) = mean(varianceInLoopMMA);

         varianceMIndividual(k) = mean(varianceInLoopIndividual);
         varianceMEqual1(k) = averagingVariance(equalWeights1, psiMatrixTrue);
         varianceMEqual2(k) = averagingVariance(equalWeights2, psiMatrixTrue);
         varianceMEqualAll(k) = averagingVariance(equalWeightsAll, psiMatrixTrue);    

    end
     
       varianceSample1(j) = mean(varianceMSample1./varianceMIndividual);
     varianceSample2(j) = mean(varianceMSample2./varianceMIndividual);
     varianceSampleAll(j) = mean(varianceMSampleAll./varianceMIndividual);
     varianceStein1(j) =  mean(varianceMStein1./varianceMIndividual);
     varianceStein2(j) =  mean(varianceMStein2./varianceMIndividual);
     varianceStein(j) =  mean(varianceMStein./varianceMIndividual);
     varianceAIC1(j) = mean(varianceMAIC1./varianceMIndividual);
     varianceAIC2(j) = mean(varianceMAIC2./varianceMIndividual);
     varianceAIC(j) = mean(varianceMAIC./varianceMIndividual);
     varianceMMA1(j) = mean(varianceMMMA1./varianceMIndividual);
     varianceMMA2(j) = mean(varianceMMMA2./varianceMIndividual);
     varianceMMA(j) = mean(varianceMMMA./varianceMIndividual);

     varianceIndividual(j) = mean(varianceMIndividual./varianceMIndividual);
     varianceEqual1(j) =  mean(varianceMEqual1./varianceMIndividual);
     varianceEqual2(j) =  mean(varianceMEqual2./varianceMIndividual);
     varianceEqualAll(j) =  mean(varianceMEqualAll); 
     
end

%  Plot 2
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
%     ylim([0q, 2])
    xlabel('\lambda_1')

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
    xlabel('\lambda_1')

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
    xlabel('\lambda_1')
    ylim([0, 2])
% 
suptitle('Averaging Estimator, \mu(\theta_1) = \lambda_1, ratio of MSE to individual estimator')
% suptitle('Averaging Estiator, \mu(\theta_1) = E(y_{T+1}|y_T=1=-1, x=1), ratio of MSE to individual estimator')
%  suptitle('Averaging Estimator, \mu(\theta_1) = {\beta_1}/{1-\lambda_1}, ratio of MSE to individual estimator')


%% 3. Coverage for plug-in CI
