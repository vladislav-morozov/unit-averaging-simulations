
% Draw coefficients and coverages
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

 
 % Create vectors to house coverages
    % Individual
 coverageIndividual = zeros(endOfRange, 1);
    % Plug-in optimized
 coverageSample1 = zeros(endOfRange, 1);
 coverageSample2 = zeros(endOfRange, 1);
 coverageSampleAll = zeros(endOfRange, 1);
    % Population optimal
 coverageP1 = zeros(endOfRange, 1);
 coverageP2 = zeros(endOfRange, 1);
 coveragePAll = zeros(endOfRange, 1);
    % Stein-type
 coverageStein1 = zeros(endOfRange, 1);
 coverageStein2 = zeros(endOfRange, 1);
 coverageStein = zeros(endOfRange, 1);
    % Equal
 coverageEqual1 = zeros(endOfRange, 1);
 coverageEqual2 = zeros(endOfRange, 1);
 coverageEqualAll = zeros(endOfRange, 1);
    % AIC/BIC
 coverageAIC1 = zeros(endOfRange, 1);
 coverageAIC2 = zeros(endOfRange, 1);
 coverageAIC = zeros(endOfRange, 1);
    % MMA
 coverageMMA1 = zeros(endOfRange, 1);
 coverageMMA2 = zeros(endOfRange, 1);
 coverageMMA = zeros(endOfRange, 1);

 
 % Length
 lengthIndivual =   zeros(endOfRange, 1);
 
 lengthSample1  =   zeros(endOfRange, 1);
 lengthSample2 =   zeros(endOfRange, 1);
 lengthSampleAll  =   zeros(endOfRange, 1);
 lengthEqual1  =   zeros(endOfRange, 1);
 lengthEqual2 =   zeros(endOfRange, 1);
 lengthEqualAll  =   zeros(endOfRange, 1);
 
% Main loop: loop through parameter vector, draw multiple samples
for j=1:endOfRange
   
   % Change coordinate
   betaLoop =    dynamicChangeCoordinate1(lambdaRange(j), betaOriginal);
    
   % Gradient describes the parameter of interest
   D = [1; 0]; % first coefficient
%    muLoop = betaLoop(1, 1);
%     D = [betaLoop(2, 1)/(1-betaLoop(1, 1))^2; 1/(1-betaLoop(1, 1))]; % LR
%     D = [2*betaLoop(2,1); 2*betaLoop(1, 1)];  % norm of beta
 
   etaTrue = betaLoop - repmat(meanCoef, 1, N); % obtain deviations from the mean
   V = linearDynamicTrueCovariance(betaLoop(1,:), betaLoop(2,:),sigmaSq, N,design);    % compute population-fixed variances
   [psiMatrixTrue, ~, ~] = averagingPopulationPsi(etaTrue, V, D); % recompute population Psi
%    psiMatrixTrue = averagingPopulationPsiOLD(etaTrue, V, D);
   [~, psiStein1] = averagingSteinWeights(etaTrue, compValue1, V, D); % create true matrix for Stein
   [~, psiStein2] = averagingSteinWeights(etaTrue, compValue2, V, D); % create true matrix for Stein
   [~, psiStein] = averagingSteinWeights(etaTrue, N, V, D); % create true matrix for Stein

   % Recreate coverage vectors
   coverageInLoopIndividual = zeros(numReplications, 1);
   coverageInLoopSample1 = zeros(numReplications, 1);
   coverageInLoopSample2 = zeros(numReplications, 1);
   coverageInLoopSampleAll = zeros(numReplications, 1);
   coverageInLoopStein1 = zeros(numReplications, 1);
   coverageInLoopAIC1 = zeros(numReplications, 1);
   coverageInLoopAIC2 = zeros(numReplications, 1);
   coverageInLoopAIC = zeros(numReplications, 1);
   coverageInLoopMMA1 = zeros(numReplications, 1);
   coverageInLoopMMA2 = zeros(numReplications, 1);
   coverageInLoopMMA = zeros(numReplications, 1);
   coverageInLoopOptimalP1 = zeros(numReplications, 1);
   coverageInLoopOptimalP2 = zeros(numReplications, 1);
   coverageInLoopOptimalPAll = zeros(numReplications, 1);
   coverageInLoopIndividual = zeros(numReplications, 1);
   coverageInLoopEqual1 =   zeros(numReplications, 1);
   coverageInLoopEqual2 =   zeros(numReplications, 1);
   coverageInLoopEqualAll =   zeros(numReplications, 1);
   
   lengthInIndividual =   zeros(numReplications, 1);
   
   lengthInSample1=   zeros(numReplications, 1);
   lengthInSample2=   zeros(numReplications, 1);
   lengthInSampleAll=   zeros(numReplications, 1);
   lengthInEqual1=   zeros(numReplications, 1);
   lengthInEqual2=   zeros(numReplications, 1);
   lengthInEqualAll=   zeros(numReplications, 1);
   
   % Create new vector to fill with MC replications
    
   betaScaled = meanCoef+etaTrue/sqrt(T);
%    betaScaled = dynamicChangeCoordinate1(lambdaRange(j), betaScaled);
        
    parfor i=1:numReplications
        
        
        [y, x, u] = linearDynamicSimulateData(betaScaled, sigmaSq, T, design, i); % draw data
        betaHat = linearStaticEstimators(y, x);% estimate coefficients
        etaEstimated = sqrt(T)*(betaHat - repmat(mean(betaHat,2),1, N));
        Vest = linearDynamicVarianceEstimator(y, x, betaHat); % estimate variance
        Dest= D; % plug into the gradient
        muLoop = betaScaled(1, 1);
        muHat = betaHat(1,:);
%          Dest = [2*betaHat(1,1); 2*betaHat(2,1)];
%         Dest = [betaHat(2, 1)/(1-betaHat(1, 1))^2; 1/(1-betaHat(1, 1))];
    
        psiMatrixSample = averagingSamplePsi(betaHat*sqrt(T), Vest, Dest);% estimate sample psi
        [weightsStein1, psiSteinSample1] = averagingSteinWeights(betaHat*sqrt(T), compValue1, Vest, Dest);
        [weightsStein2, psiSteinSample2] = averagingSteinWeights(betaHat*sqrt(T), compValue2, Vest, Dest);
        [weightsStein, psiSteinSample] = averagingSteinWeights(betaHat*sqrt(T), N, Vest, Dest);

        % Individual
        wInd = zeros(N, 1);
        wInd(1) = 1;
        [bSI, vSI] = averagingNaiveCIbounds(betaHat, wInd, Dest, Vest, T);
        coverageInLoopIndividual(i) = muLoop>= betaHat(1,1)-bSI-vSI...
                            && muLoop<= betaHat(1,1)-bSI+vSI;
        lengthInIndividual(i) = 2*vSI;
        
        % Sample block
        weightsSample1 = averagingOptimalWeights(psiMatrixSample(1:compValue1,1:compValue1), nW);
        muHatSample1 = muHat(1:compValue1)*weightsSample1;
        [bS1, vS1] = averagingNaiveCIbounds(betaHat, weightsSample1, Dest, Vest, T);
 
        coverageInLoopSample1(i) = muLoop>= muHatSample1-bS1-vS1...
                            && muLoop<= muHatSample1-bS1+vS1;
        lengthInSample1(i) = 2*vS1;
        
        weightsSample2 = averagingOptimalWeights(psiMatrixSample(1:compValue2,1:compValue2), nW);
        muHatSample2 = muHat(1:compValue2)*weightsSample2;
        [bS2, vS2] = averagingNaiveCIbounds(betaHat, weightsSample2, Dest, Vest, T);
 
        coverageInLoopSample2(i) = muLoop>= muHatSample2-bS2-vS2...
                            && muLoop<= muHatSample2-bS2+vS2;
        lengthInSample2(i) = 2*vS2;
         
        weightsSampleAll = averagingOptimalWeights(psiMatrixSample(1:N,1:N), nW);
        muHatSampleAll = muHat(1:N)*weightsSampleAll;
        [bSAll, vSAll] = averagingNaiveCIbounds(betaHat, weightsSampleAll, Dest, Vest, T);
 
        coverageInLoopSampleAll(i) = muLoop>= muHatSampleAll-bSAll-vSAll...
                            && muLoop<= muHatSampleAll-bSAll+vSAll;
        lengthInSampleAll(i) = 2*vSAll;
        
        
        
        % Equal weights
        
         muHatEqual1 = muHat*equalWeights1;
        [bE1, vE1] = averagingNaiveCIbounds(betaHat, equalWeights1, Dest, Vest, T);
 
        coverageInLoopEqual1(i) = muLoop>= muHatEqual1-bE1-vE1...
                            && muLoop<= muHatEqual1-bE1+vE1;
        lengthInEqual1(i) = 2*vE1;
        
                 muHatEqual2 = muHat*equalWeights2;
        [bE2, vE2] = averagingNaiveCIbounds(betaHat, equalWeights2, Dest, Vest, T);
 
        coverageInLoopEqual2(i) = muLoop>= muHatEqual1-bE2-vE2...
                            && muLoop<= muHatEqual2-bE2+vE2;
        lengthInEqual2(i) = 2*vE2;
        
                 muHatEqualAll = muHat*equalWeightsAll;
        [bEAll, vEAll] = averagingNaiveCIbounds(betaHat, equalWeightsAll, Dest, Vest, T);
        coverageInLoopEqualAll(i) = muLoop>= muHatEqualAll-bEAll-vEAll...
                            && muLoop<= muHatEqualAll-bEAll+vEAll;
        lengthInEqualAll(i) = 2*vEAll;
        [j, i]
    end
     coverageIndividual(j)  = mean(coverageInLoopIndividual);
     coverageSample1(j) = mean(coverageInLoopSample1);
     coverageSample2(j) = mean(coverageInLoopSample2);
     coverageSampleAll(j) = mean(coverageInLoopSampleAll);
     coverageStein1(j) =  mean(coverageInLoopStein1);
     coverageStein2(j) =  mean(coverageInLoopStein1);
     coverageStein(j) =  mean(coverageInLoopStein1);
     coverageAIC1(j) = mean(coverageInLoopAIC1);
     coverageAIC2(j) = mean(coverageInLoopAIC2);
     coverageAIC(j) = mean(coverageInLoopAIC);
     coverageMMA1(j) = mean(coverageInLoopMMA1);
     coverageMMA2(j) = mean(coverageInLoopMMA2);
     coverageMMA(j) = mean(coverageInLoopMMA);

     coverageIndividual(j) = mean(coverageInLoopIndividual);
     coverageEqual1(j) =  mean(coverageInLoopEqual1);
     coverageEqual2(j) = mean(coverageInLoopEqual2);
     coverageEqualAll(j) = mean(coverageInLoopEqualAll);    
     
     
     
     lengthIndivual(j) = mean(lengthInIndividual);
     
     lengthSample1(j) = mean(lengthInSample1);
     lengthSample2(j)  = mean(lengthInSample2);
     lengthSampleAll(j)  = mean(lengthInSampleAll);
     
     lengthEqual1(j) = mean(lengthInEqual1);
     lengthEqual2(j)  = mean(lengthInEqual2);
     lengthEqualAll(j)  = mean(lengthInEqualAll);
end

%

%%
figure

subplot(1,3,1)
    
    plot(lambdaRange,coverageIndividual) % optimal
    hold on
    plot(lambdaRange,coverageSample1,'c','LineWidth',2) % optimal
    plot(lambdaRange,coverageEqual1,'--')

    legend('Individual','Plug-In', 'Equal','Stein-like', 'AIC','MMA', 'Location', 'northeast')
    title('N=5')
    xlim([min(lambdaRange), max(lambdaRange)])
    ylim([0.7, 1])
%     xlabel('\eta_1')
xlabel('\eta_1')
    
subplot(1,3,2)
    plot(lambdaRange,coverageIndividual) % optimal
    hold on
    plot(lambdaRange,coverageSample2,'c','LineWidth',2) % optimal
    plot(lambdaRange,coverageEqual2,'--')

    legend('Individual','Plug-In', 'Equal','Stein-like', 'AIC','MMA', 'Location', 'northeast')
    title('N=25')
    xlim([min(lambdaRange), max(lambdaRange)])
    ylim([0.7, 1])
    xlabel('\eta_1')
%     xlabel('\lambda_1')
subplot(1,3,3)
    plot(lambdaRange,coverageIndividual) % optimal
    hold on
    plot(lambdaRange,coverageSampleAll,'c','LineWidth',2) % optimal
    plot(lambdaRange,coverageEqualAll,'--')

    legend('Individual','Plug-In', 'Equal','Stein-like', 'AIC','MMA', 'Location', 'northeast')
    title('N=100')
    xlim([min(lambdaRange), max(lambdaRange)])
    ylim([0.7, 1])
    xlabel('\eta_1')
%     xlabel('\lambda_1')
suptitle('Coverage probabilities, 1-step plug-in CI, Nominal 95%')

%  Average length
figure

subplot(1,3,1)
    plot(lambdaRange,lengthIndivual,'LineWidth',2) % optimal
    hold on
    plot(lambdaRange,lengthEqual1,'--')
    plot(lambdaRange,lengthSample1,'--')
    legend('Individual', 'Plug-In', 'Equal','Stein-like', 'AIC','MMA', 'Location', 'northeast')
    title('N=5')
    xlim([min(lambdaRange), max(lambdaRange)])
    ylim([0, 0.7])
    xlabel('\eta_1')
% xlabel('\lambda_1')
    
subplot(1,3,2)
    plot(lambdaRange,lengthIndivual,'LineWidth',2) % optimal
    hold on
    plot(lambdaRange,lengthEqual2,'--')
    plot(lambdaRange,lengthSample2,'--')
    legend('Individual', 'Plug-In', 'Equal','Stein-like', 'AIC','MMA', 'Location', 'northeast')
    title('N=25')
    xlim([min(lambdaRange), max(lambdaRange)])
    ylim([0, 0.7])
    xlabel('\eta_1')
%     xlabel('\lambda_1')
subplot(1,3,3)
    plot(lambdaRange,lengthIndivual,'LineWidth',2) % optimal
    hold on
    plot(lambdaRange,lengthEqualAll,'--')
    plot(lambdaRange,lengthSampleAll,'--')
    legend('Individual', 'Plug-In', 'Equal','Stein-like', 'AIC','MMA', 'Location', 'northeast')
    title('N=100')
    xlim([min(lambdaRange), max(lambdaRange)])
    ylim([0, 0.7])
    xlabel('\eta_1')
%     xlabel('\lambda_1')
suptitle('Length, 1-step plug-in CI, Nominal 95%')

%% Average length