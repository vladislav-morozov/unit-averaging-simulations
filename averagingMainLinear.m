
%% Set some parameters

N = 100;
T =6;
design = 1;
meanBeta = 0;

%% Create data
seedData = 14;
seedCoefficients = 14;
rng(seedData,'philox')
[y0, x, u, lambda, betaC, lambdaScaled, betaScaled, sigmaSq] =...
        linearDynamicDrawCoefficients(N, T, design, meanBeta, seedData, seedCoefficients);
eta = [lambdaScaled; betaScaled];
y = linearDynamicSimulateData(eta, y0, x, u);

%% Get individual estimators and the MG estimator
theta = linearDynamicEstimators(y, x, N, T);
thetaMG = mean(theta);
% Obtain covariances, use true non-scaled coefficients
V = linearDynamicTrueCovariance(   lambda, betaC, sigmaSq, N);



%% Plot asymptotic MSE
% Note: use design = 1, scaling sigmaSq by 10 times produces a more compact
% plot

% Suppose we wish to estimate the beta coefficient
lambda1 = 0.25:0.01:0.75; % range of values for the first unit
D = [1; 0]; % that we only care about the second parameter

% Compute variances quickly using given weights and varying the parameter
% of the first uni. Use all, first 5 and first 25 units
varianceWeighted = @(w, l1) averagingVariance(w,...
         averagingPopulationPsi( [[l1 lambda(2:end)]-0.5; betaC(1:end)],V,D));
varianceWeightedStein = @(w, l1) averagingVariance(w,...
         averagingPopulationPsi( [[l1]-0.5; betaC(1:1)],V,D));
varianceWeighted5 = @(w, l1) averagingVariance(w,...
         averagingPopulationPsi( [[l1 lambda(2:5)]-0.5; betaC(1:5)],V,D));
varianceWeighted25 = @(w, l1) averagingVariance(w,...
         averagingPopulationPsi( [[l1 lambda(2:25)]-0.5; betaC(1:25)],V,D));
     
% Compute optimal weights given coefficient for unit 1
nW = 0; % allow/disallow negative weights to occur

optimalWeights = @(l1) averagingOptimalWeights( ...
                averagingPopulationPsi( [[l1 lambda(2:end)]-0.5; betaC(1:end)],V,D), nW);
optimalWeightsStein = @(l1) averagingOptimalWeights( ...
                averagingPopulationPsi( [l1; betaC(1:1)],V,D), nW);            
optimalWeights5 = @(l1) averagingOptimalWeights( ...
                averagingPopulationPsi( [[l1 lambda(2:5)]-0.5; betaC(1:5)],V,D), nW);
optimalWeights25 = @(l1) averagingOptimalWeights( ...
                averagingPopulationPsi( [[l1 lambda(2:25)]-0.5; betaC(1:25)],V,D), nW);
           

                                    
% Create equal weights            
equalWeights5 = zeros(N+1,1);
equalWeights5(1:6) = ones(6,1)/6;
equalWeightsAll = ones(N+1,1)/(N+1);

% Create plot
figure
plot(lambda1, linearDynamicVarianceLambda1(lambda1, betaC(1), sigmaSq(1))); % AMSE of individual
hold on
plot(lambda1, (lambda1-0.5).^2);  % AMSE of MG
plot(lambda1, arrayfun(@(x) varianceWeighted(equalWeights5, x), lambda1)) % equal weights, N=5
plot(lambda1, arrayfun(@(x) varianceWeighted(equalWeightsAll, x), lambda1)) % equal weights, N=25

% Optimal
plot(lambda1, arrayfun(@(x) varianceWeightedStein(optimalWeightsStein(x), x), lambda1),'--g','LineWidth',2) % optimal
plot(lambda1, arrayfun(@(x) varianceWeighted5(optimalWeights5(x), x), lambda1),'c','LineWidth',2) % optimal
plot(lambda1, arrayfun(@(x) varianceWeighted25(optimalWeights25(x), x), lambda1),'k','LineWidth',2) % optimal
plot(lambda1, arrayfun(@(x) varianceWeighted(optimalWeights(x), x), lambda1),'m','LineWidth',2) % optimal



legend('Individual', 'Mean group', 'Equal weights, N=5', 'Equal weights, N=25', ...
    'Optimal Weights, Stein-type', 'Optimal Weights, N=5','Optimal Weights, N=25','Optimal Weights, N=100')
xlabel('\lambda')
% xlabel('\eta_1, Difference of coefficient of unit 1 from population mean ')
xlim([min(lambda1), max(lambda1)])
ylabel('MSE')
 ylim([-.01, 0.15])


%% Now comparing things in-sample

% Note: use design = 1, scaling sigmaSq by 10 times produces a more compact
% plot

% Suppose we wish to estimate the beta coefficient
lambda1 = 0.25:0.01:0.75; % range of values for the first unit
D = [1;0]; % that we only care about the second parameter

% Resimulate data with different parameter value, 
resimulateUnit1y =  @(l1) linearDynamicSimulateData([l1/sqrt(T); betaScaled(1)], y0(1), x(:,1), u(:,1));
recomputeUnit1Theta = @(l1) linearDynamicEstimators(resimulateUnit1y(l1), x(:,1), 1, T);
varianceEstimator = @(l1) cat(3, linearDynamicVarianceEstimator(resimulateUnit1y(l1),...
        x(:,1),recomputeUnit1Theta(l1)), V(:,:,2:end));

     
% Compute optimal weights given coefficient for unit 1
nW = 0; % allow/disallow negative weights to occur

optimalSampleWeights = @(l1) averagingOptimalWeights( ...
                averagingSamplePsi([recomputeUnit1Theta(l1), theta(:,2:end)],...
                varianceEstimator(l1),D), nW);
optimalSampleWeightsStein = @(l1) averagingOptimalWeights( ...
                averagingSamplePsi( [l1; betaC(1:1)],V,D), nW); 
optimalSampleWeights5 = @(l1) averagingOptimalWeights( ...
                averagingSamplePsi([recomputeUnit1Theta(l1), theta(:,2:5)],...
                varianceEstimator(l1),D), nW);
optimalSampleWeights25 = @(l1) averagingOptimalWeights( ...
                averagingSamplePsi([recomputeUnit1Theta(l1), theta(:,2:25)],...
                varianceEstimator(l1),D), nW);

% Create equal weights            
equalWeights5 = zeros(N+1,1);
equalWeights5(1:6) = ones(6,1)/6;
equalWeightsAll = ones(N+1,1)/(N+1);

% Create plot
figure
plot(lambda1, linearDynamicVarianceLambda1(lambda1, betaC(1), sigmaSq(1))); % AMSE of individual
hold on
plot(lambda1, (lambda1-0.5).^2);  % AMSE of MG
plot(lambda1, arrayfun(@(x) varianceWeighted(equalWeights5, x), lambda1)) % equal weights, N=5
plot(lambda1, arrayfun(@(x) varianceWeighted(equalWeightsAll, x), lambda1)) % equal weights, N=25

% Optimal
plot(lambda1, arrayfun(@(x) varianceWeightedStein(optimalSampleWeightsStein(x), x), lambda1),'--g','LineWidth',2) % optimal
plot(lambda1, arrayfun(@(x) varianceWeighted5(optimalSampleWeights5(x), x), lambda1),'c','LineWidth',2) % optimal
plot(lambda1, arrayfun(@(x) varianceWeighted25(optimalSampleWeights25(x), x), lambda1),'k','LineWidth',2) % optimal
plot(lambda1, arrayfun(@(x) varianceWeighted(optimalSampleWeights(x), x), lambda1),'m','LineWidth',2) % optimal



legend('Individual', 'Mean group', 'Equal weights, N=5', 'Equal weights, N=25', ...
    'Optimal Plug-in Weights, Stein-type', 'Optimal Weights, N=5',...
    'Optimal Plug-in Weights, N=25','Optimal Plug-in Weights, N=100')
xlabel('\lambda')
% xlabel('\eta_1, Difference of coefficient of unit 1 from population mean ')
xlim([min(lambda1), max(lambda1)])
ylabel('MSE')
hold off
 ylim([-.01, 0.15])
 
 
 %% Parfor
 numReplications  = 100; % number of replications
 
 lambdaTried = 0.25:0.01:0.75;
 nLambda = length(lambdaTried);
 
 nW = 0;
 varianceIndividual = zeros(nLambda, 1);
 varianceSample5 = zeros(nLambda, 1);
 varianceSample25 = zeros(nLambda, 1);
 varianceSampleAll = zeros(nLambda, 1);
 varianceP5 = zeros(nLambda, 1);
 varianceP25 = zeros(nLambda, 1);
 variancePAll = zeros(nLambda, 1);
 varianceEqual5 = zeros(nLambda, 1);
 varianceEqual25 = zeros(nLambda, 1);
 varianceEqualAll = zeros(nLambda, 1);
 
 % Create equal weights
 equalWeights5 = zeros(N+1,1);
 equalWeights5(2:6) = ones(5,1)/5;
 equalWeights25 = zeros(N+1,1);
 equalWeights25(2:26) = ones(25,1)/25;
 equalWeightsAll = ones(N+1,1)/(N+1);
 equalWeightsAll(1) =0;
 
 % Draw coefficients
 meanBeta = 2;
 [~, ~, ~, lambda, betaC, lambdaScaled, betaScaled, sigmaSq] =...
            linearDynamicDrawCoefficients(N, T, design, meanBeta, 14, 16);
        % 16
   
 % Compute for equal weights

 for j=1:nLambda
     l1 = lambdaTried(j);
     varianceInLoopOptimal5 = zeros(numReplications, 1);
     varianceInLoopOptimal25 = zeros(numReplications, 1);
     varianceInLoopOptimalAll = zeros(numReplications, 1);
     varianceInLoopOptimalP5 = zeros(numReplications, 1);
     varianceInLoopOptimalP25 = zeros(numReplications, 1);
     varianceInLoopOptimalPAll = zeros(numReplications, 1);
     varianceInLoopIndividual = zeros(numReplications, 1);
     lambda(1) = l1;
     
     etaTrue = [lambda-0.5; betaC-meanBeta];
     % Obtain covariances, use true non-scaled coefficients
     V = linearDynamicTrueCovariance(   lambda, betaC, sigmaSq, N);
     % Recompute population matrices which do not depend on data   
     psiMatrixTrue5 = averagingPopulationPsi(etaTrue(:,1:6) , V, D);
     psiMatrixTrue25 = averagingPopulationPsi(etaTrue(:,1:26) , V, D);
     psiMatrixTrueAll = averagingPopulationPsi(etaTrue , V, D);
     
     D = [1;0]; 
     %D = [2*lambda(1); 2*betaC(1)];
     parfor i=1:numReplications
         % Draw new data
        [y0, x, u, ~, ~, ~, ~, ~] =...
            linearDynamicDrawCoefficients(N, T, design, meanBeta, i, i);
        
        
%         betaC(1) = 1.28;
%         sigmaSq(1) = 0.1;
        coefScaled = [lambdaScaled; betaScaled];
        y = linearDynamicSimulateData(coefScaled, y0, x, u);
        theta = linearDynamicEstimators(y, x, N, T);
        thetaMG = mean(theta);
        Vestimated = linearDynamicVarianceEstimator(y, x, theta);
        
        
        % Compute sample
        psiMatrixSample5 = averagingSamplePsi(theta(:, 1:6), Vestimated, D);   
        psiMatrixSample25 = averagingSamplePsi(theta(:, 1:26), Vestimated, D);
        psiMatrixSampleAll = averagingSamplePsi(theta, Vestimated, D);
     
        
        varianceInLoopOptimalP5(i) = averagingVariance(...
                averagingOptimalWeights(psiMatrixTrue5, nW), psiMatrixTrue5); 
        varianceInLoopOptimalP25(i) = averagingVariance(...
                averagingOptimalWeights(psiMatrixTrue25, nW), psiMatrixTrue25); 
        varianceInLoopOptimalPAll(i) = averagingVariance(...
                averagingOptimalWeights(psiMatrixTrueAll, nW), psiMatrixTrueAll);
        varianceInLoopOptimal5(i) = averagingVariance(...
                averagingOptimalWeights(psiMatrixSample5, nW), psiMatrixTrue5); 
        varianceInLoopOptimal25(i) = averagingVariance(...
                averagingOptimalWeights(psiMatrixSample25, nW), psiMatrixTrue25); 
        varianceInLoopOptimalAll(i) = averagingVariance(...
                averagingOptimalWeights(psiMatrixSampleAll, nW), psiMatrixTrueAll); 
        varianceInLoopIndividual(i) = D'*V(:, :, 1)*D;    
            
        
        [j, i]
     end
     varianceSample5(j) = mean(varianceInLoopOptimal5);
     varianceSample25(j) = mean(varianceInLoopOptimal25);
     varianceSampleAll(j) = mean(varianceInLoopOptimalAll);
     varianceP5(j) = mean(varianceInLoopOptimalP5);
     varianceP25(j) = mean(varianceInLoopOptimalP25);
     variancePAll(j) = mean(varianceInLoopOptimalPAll);
     varianceIndividual(j) = mean(varianceInLoopIndividual);
     varianceEqual5(j) = averagingVariance(equalWeights5, psiMatrixTrueAll);
     varianceEqual25(j) = averagingVariance(equalWeights25, psiMatrixTrueAll);
     varianceEqualAll(j) = averagingVariance(equalWeightsAll, psiMatrixTrueAll);
 end
 
    
 figure
%  plot(lambdaTried,varianceIndividual) % individual
 
%  plot(lambdaTried,(lambdaTried-0.5).^2./varianceIndividual')  % MG
 hold on
 plot(lambdaTried,varianceSample5./varianceIndividual)
  plot(lambdaTried,varianceSample25./varianceIndividual)
   plot(lambdaTried,varianceSampleAll./varianceIndividual)
   
 plot(lambdaTried,varianceEqual5./varianceIndividual,'--')
  plot(lambdaTried,varianceEqual25./varianceIndividual,'--')
   plot(lambdaTried,varianceEqualAll./varianceIndividual,'--')
      
%     plot(lambdaTried,varianceP5)
%   plot(lambdaTried,varianceP25)
%    plot(lambdaTried,variancePAll)
   
 legend( 'N=5', 'N=25', 'N=100', 'Equal 5', 'Equal 25', 'Equal 100')
