
%% Set some parameters

N = 100;
T =10;
design = 1;
meanBeta = 0;

%% Create data
seed = 14;
rng(seed,'philox')
[y0, x, u, lambda, beta, lambdaScaled, betaScaled, sigmaSq] =...
        linearDynamicDrawCoefficients(N, T, design, meanBeta, seed);
eta = [lambdaScaled; betaScaled];
y = linearDynamicSimulateData(eta, y0, x, u);

%% Get individual estimators and the MG estimator
theta = linearDynamicEstimators(y, x, N, T);
thetaMG = mean(theta);
% Obtain covariances, use true non-scaled coefficients
V = linearDynamicTrueCovariance(   lambda, beta, sigmaSq, N);



%% Plot asymptotic MSE
% Note: use design = 1, scaling sigmaSq by 10 times produces a more compact
% plot

% Suppose we wish to estimate the beta coefficient
beta1 = -5:0.01:5; % range of values for the first unit
D = [0;1]; % that we only care about the second parameter

% Compute variances quickly using given weights and varying the parameter
% of the first uni. Use all, first 5 and first 25 units
varianceWeighted = @(w, l1) averagingVariance(w,...
         averagingPopulationPsi([ lambda(1:end);[l1, beta(2:end)]],V,D));
varianceWeighted5 = @(w, l1) averagingVariance(w,...
         averagingPopulationPsi([ lambda(1:5);[l1, beta(2:5)]],V,D));
varianceWeighted25 = @(w, l1) averagingVariance(w,...
         averagingPopulationPsi([ lambda(1:25);[l1, beta(2:25)]],V,D));
     
% Compute optimal weights given coefficient for unit 1
nW = 0; % allow/disallow negative weights to occur

optimalWeights = @(l1) averagingOptimalWeights( ...
                averagingPopulationPsi([ lambda(1:end);[l1, beta(2:end)]],V,D), nW);
optimalWeights5 = @(l1) averagingOptimalWeights( ...
                averagingPopulationPsi([ lambda(1:5);[l1, beta(2:5)]],V,D), nW);
optimalWeights25 = @(l1) averagingOptimalWeights( ...
                averagingPopulationPsi([ lambda(1:25);[l1, beta(2:25)]],V,D), nW);
           
optimalSampleWeights = @(l1) averagingOptimalWeights( ...
                averagingPopulationPsi([ lambda(1:end);[l1, beta(2:end)]],V,D));
optimalSampleWeights5 = @(l1) averagingOptimalWeights( ...
                averagingPopulationPsi([ lambda(1:5);[l1, beta(2:5)]],V,D));
optimalSampleWeights25 = @(l1) averagingOptimalWeights( ...
                averagingPopulationPsi([ lambda(1:25);[l1, beta(2:25)]],V,D));
                                    
% Create equal weights            
equalWeights5 = zeros(N+1,1);
equalWeights5(1:6) = ones(6,1)/6;
equalWeightsAll = ones(N+1,1)/(N+1);

% Create plot
figure
plot(beta1, linearDynamicVarianceLambda2(lambda(1), beta1, sigmaSq(1))); % AMSE of individual
hold on
plot(beta1, beta1.^2);  % AMSE of MG
plot(beta1, arrayfun(@(x) varianceWeighted(equalWeights5, x), beta1)) % equal weights, N=5
plot(beta1, arrayfun(@(x) varianceWeighted(equalWeightsAll, x), beta1)) % equal weights, N=25

% Optimal
plot(beta1, arrayfun(@(x) varianceWeighted5(optimalWeights5(x), x), beta1),'c','LineWidth',2) % optimal
plot(beta1, arrayfun(@(x) varianceWeighted25(optimalWeights25(x), x), beta1),'k','LineWidth',2) % optimal
plot(beta1, arrayfun(@(x) varianceWeighted(optimalWeights(x), x), beta1),'m','LineWidth',2) % optimal



legend('Individual', 'Mean group', 'Equal weights, N=5', 'Equal weights, N=25', ...
    'Optimal Weights, N=5','Optimal Weights, N=25','Optimal Weights, N=100')
xlabel('\eta')
% xlabel('\eta_1, Difference of coefficient of unit 1 from population mean ')
xlim([min(beta1), max(beta1)])
ylabel('MSE')
 ylim([-.1, 1])


%% Now comparing things in-sample

% Note: use design = 1, scaling sigmaSq by 10 times produces a more compact
% plot

% Suppose we wish to estimate the beta coefficient
beta1 = -5:0.01:5; % range of values for the first unit
D = [0;1]; % that we only care about the second parameter

% Resimulate data with different parameter value, 
resimulateUnit1y =  @(l1) linearDynamicSimulateData([lambdaScaled(1); l1/sqrt(T)], y0(1), x(:,1), u(:,1));
recomputeUnit1Theta = @(l1) linearDynamicEstimators(resimulateUnit1y(l1), x(:,1), 1, T);
varianceEstimator = @(l1) cat(3, linearDynamicVarianceEstimator(resimulateUnit1y(l1),...
        x(:,1),recomputeUnit1Theta(l1)), V(:,:,2:end));

     
% Compute optimal weights given coefficient for unit 1
nW = 1; % allow/disallow negative weights to occur

optimalSampleWeights = @(l1) averagingOptimalWeights( ...
                averagingSamplePsi([recomputeUnit1Theta(l1), theta(:,2:end)],...
                varianceEstimator(l1),D), nW);
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
plot(beta1, linearDynamicVarianceLambda2(lambda(1), beta1, sigmaSq(1))); % AMSE of individual
hold on
plot(beta1, beta1.^2);  % AMSE of MG
plot(beta1, arrayfun(@(x) varianceWeighted(equalWeights5, x), beta1)) % equal weights, N=5
plot(beta1, arrayfun(@(x) varianceWeighted(equalWeightsAll, x), beta1)) % equal weights, N=25

% Optimal
plot(beta1, arrayfun(@(x) varianceWeighted5(optimalSampleWeights5(x), x), beta1),'c','LineWidth',2) % optimal
plot(beta1, arrayfun(@(x) varianceWeighted25(optimalSampleWeights25(x), x), beta1),'k','LineWidth',2) % optimal
plot(beta1, arrayfun(@(x) varianceWeighted(optimalSampleWeights(x), x), beta1),'m','LineWidth',2) % optimal



legend('Individual', 'Mean group', 'Equal weights, N=5', 'Equal weights, N=25', ...
    'Optimal Weights, N=5','Optimal Plug-in Weights, N=25','Optimal Plug-in Weights, N=100')
xlabel('\eta')
% xlabel('\eta_1, Difference of coefficient of unit 1 from population mean ')
xlim([min(beta1), max(beta1)])
ylabel('MSE')
 ylim([-.1, 1])
