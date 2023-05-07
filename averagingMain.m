%% MSE Illustration

% Set some parameters

% Create data
rng(14,'philox')
[y, x, lambda, beta] = linearDynamicSimulateData(N, T, design);

% Get individual estimators and the MG estimator
theta = linearDynamicEstimators(y, x, N, T);
thetaMG = mean(theta);