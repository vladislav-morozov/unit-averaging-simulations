

%% Define parameters

compValue1 = 5;
compValue2 = 25;

seedData = 16;
seedCoefficients =15;
meanBeta = [0;0];
varianceBeta = 1;
correlationBeta = 0.6;
design =1;

endOfRange = 5;
N= 100;
T = 10;

numReplications = 10;
nW = 0;
%% Create objects for later use
 % Create equal weights
 equalWeights1 = zeros(N,1);
 equalWeights1(1:compValue1) = ones(compValue1,1)/5;
 equalWeights2 = zeros(N,1);
 equalWeights2(1:compValue2) = ones(compValue2,1)/compValue2;
 equalWeightsAll = ones(N,1)/(N);
%  equalWeightsAll(1) =0;
 
 varianceIndividual = zeros(endOfRange, 1);
 varianceSample1 = zeros(endOfRange, 1);
 varianceSample2 = zeros(endOfRange, 1);
 varianceSampleAll = zeros(endOfRange, 1);
 varianceP1 = zeros(endOfRange, 1);
 varianceP2 = zeros(endOfRange, 1);
 variancePAll = zeros(endOfRange, 1);
 varianceEqual1 = zeros(endOfRange, 1);
 varianceEqual2 = zeros(endOfRange, 1);
 varianceEqualAll = zeros(endOfRange, 1);


% Main loop

[betaOriginal, betaScaled, sigmaSq] = linearStaticDrawCoefficients(N, T, meanBeta, varianceBeta,...
            correlationBeta, seedCoefficients);

        
beta1range = -2.5:0.02:2.5;
endOfRange = length(beta1range);
%% 
for j=1:endOfRange
   
   b1 = beta1range(j);
   betaLoop =    staticChangeCoordinate1(b1, betaScaled);
   %  staticChangeCoordinate2
   %  staticChangeScale
   % select new value for parameter
   D = [1; 0    ]; % recompute the gradient D with new values
   etaTrue = betaLoop - repmat(meanBeta, 1, N); % obtain deviations from the mean
   V = linearStaticTrueCovariance(sigmaSq, N);    % compute population-fixed variances
   psiMatrixTrue = averagingPopulationPsiOLD(etaTrue , V, D); % recompute population Psi
   
   % Recreate variance vectors
   varianceInLoopSample1 = zeros(numReplications, 1);
   varianceInLoopSample2 = zeros(numReplications, 1);
   varianceInLoopSampleAll = zeros(numReplications, 1);
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
        psiMatrixSample = averagingSamplePsi(betaHat, Vest, Dest);% estimate sample psi
%         
%         varianceInLoopOptimalP1(i) = averagingVariance(...
%                 averagingOptimalWeights(psiMatrixTrue(2:compValue1+1,2:compValue1+1), nW),...
%                 psiMatrixTrue(2:compValue1+1,2:compValue1+1)); 
%         varianceInLoopOptimalP2(i) = averagingVariance(...
%                 averagingOptimalWeights(psiMatrixTrue(2:compValue2+1,2:compValue2+1), nW),...
%                 psiMatrixTrue(2:compValue2+1,2:compValue2+1)); 
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
        varianceInLoopIndividual(i) = D'*V(:, :, 1)*D;    
            
        
        [j, i]
    end
     varianceSample1(j) = mean(varianceInLoopSample1);
     varianceSample2(j) = mean(varianceInLoopSample2);
     varianceSampleAll(j) = mean(varianceInLoopSampleAll);
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

% Do the plots
figure
subplot(1,3,1)

 plot(beta1range,varianceSample1./varianceIndividual)
 hold on


 plot(beta1range,varianceEqual1./varianceIndividual,'--')
legend('Plug-In', 'Equal')
title('N=5')
xlim([min(beta1range), max(beta1range)])
 ylim([0, 4])

     subplot(1,3,2)
  plot(beta1range,varianceSample2./varianceIndividual)
hold on
  plot(beta1range,varianceEqual2./varianceIndividual,'--')
legend('Plug-In', 'Equal')
title('N=25')
xlim([min(beta1range), max(beta1range)])
 ylim([0, 4])

subplot(1,3,3)

   plot(beta1range,varianceSampleAll./varianceIndividual)
   hold on
   plot(beta1range,varianceEqualAll./varianceIndividual,'--')
legend('Plug-In', 'Equal')
title('N=100')
xlim([min(beta1range), max(beta1range)])
 ylim([0, 4])