% Simulations for "Unit Averaging for Heterogeneous Panels"
%
% Authors: Christian Brownlees, Vladislav Morozov
% Contact: Vladislav Morozov (vladislav.morozov (at) barcelonagse.eu)
%
% Last modified: 2024-03-10
%
%
% This file runs the simulations reported in the paper and the online
% appendix. 
% To replicate FINDING X, run the file IN THIS WAY
%
% 
% 

clear variables
close all
%% Simulation parameters

% whether plots should be shown or saved without showing
plotQuietly = false; 

% Select a seed to drawing coefficients
seedCoefficients = 1; 

% Number of replications
numReplications = 10;

% Saving
saveWeights = false;
saveUnrestricted = false;


%% DGP Parameters
% Data dimensions
valuesN= [10, 25, 50, 150, 250]; % values of N to compare, must be arranged in ascending order
valuesT = 60 ;

% For Dynamic model: select default mean and variance parameters for
% coefficients
lambdaMean = 0;
meanBeta = 1; % mean for beta coefficient
varianceBeta = 1; % variance for beta 
varNoiseVar = 1; % variance of heteroskedasticity
varianceX = 1; % to match the factor
 
localAsy = 1;

%% Parameters to estimate in the model
% The model is
%   y_{it} = theta_{i1} y_{it-1} + theta_{i2} x_{it} + u_{it}
%
% 
% PARAMETERS CONSIDERED:
% The methodology applies to smooth transformation mu(theta)
% This script implements 4 parameters of interest (mu):
% 1. theta_{i1} -- coefficient on the lag
% 2. theta_{i2} -- coefficient on the exogenous variable
% 3. long-run effect of a one-unit shift in x
% 4. one-step ahead out-of-sample forecast. Depends on the last value of y
%    and a value of xForecast specified in the definition
%
%
% IMPLEMENTATION NOTES: 
% -- The parameters are stored as cells in the cell array paramArray
% -- Each parameters is defined by a struct. The fields are the function
%    mu, the gradient, the name of the parameter for the plot, and the name
%    to use when saving
% -- the function mu and the gradient take both the coefficients and the
%    full data vectors y and x. 
% -- the plotting functions will use the description and name of the
%    parameters to automatically generate labels and file names of the
%    plots
% -- the script detects the length of the paramArray and loops through
%    all the parameters

% Defining parameters

% Coefficient on lag
paramArray{1}.mu = @(theta, x, y) theta(1, :); 
paramArray{1}.gradient = @(theta, x, y) [1;0]; 
paramArray{1}.plotDescr = "\mu(\theta_1) = \lambda_1";
paramArray{1}.saveName = "lambda";

% Coefficient on exogenous variable
paramArray{2}.mu = @(theta, x, y) theta(2, :);
paramArray{2}.gradient = @(theta, x, y) [0; 1]; 
paramArray{2}.plotDescr = "\mu(\theta_1) = \beta_1";
paramArray{2}.saveName = "beta";

% Long-run effect
paramArray{3}.mu = @(theta, x, y) theta(2, :)./(1-theta(1, :));
paramArray{3}.gradient = @(theta, x, y) ...
    [theta(2, 1)/(1-theta(1, 1))^2; 1/(1-theta(1, 1))]; 
paramArray{3}.plotDescr = "\mu(\theta_1) = {\beta_1}/{1-\lambda_1}";
paramArray{3}.saveName = "longRun";

% Forecast
xForecast = 1;
paramArray{4}.mu =  @(theta, x, y) ...
    theta(1, :).*y(end, 1)+theta(2, :)*xForecast; 
paramArray{4}.gradient = @(theta, x, y) [y(end, 1); 1];  
paramArray{4}.plotDescr = "\mu(\theta_1) = E(y_{T+1}|y_T, x_T=1)";
paramArray{4}.saveName = "forecast";

%% Averaging schemes

methodsArray{1}.weightFunction = ...
    @(y, covars, estCoefs, estCovars, gradientEstimateTarget, ...
      targetID, unrestrictedBooltargetId) ...
            double(1:size(estCoefs, 2) == 2)';
methodsArray{1}.shortName = 'ind';        
methodsArray{1}.longName = 'Individual';


methodsArray{2}.weightFunction = ...
    @(y, covars, estCoefs, estCovars, gradientEstimateTarget, ...
      targetID, unrestrictedBooltargetId) ...
            ones(size(estCoefs, 2), 1)/size(estCoefs, 2);
methodsArray{2}.shortName = 'mg';        
methodsArray{2}.longName = 'Mean Group';


methodsArray{3}.weightFunction = ...
    @(y, covars, estCoefs, estCovars, gradientEstimateTarget, ...
      targetID, unrestrictedBooltargetId) ...
            uaWeightsAICMMA(estCoefs, y, covars, 'aic');
methodsArray{3}.shortName = 'aic';        
methodsArray{3}.longName = 'AIC';


methodsArray{4}.weightFunction = ...
    @(y, covars, estCoefs, estCovars, gradientEstimateTarget, ...
      targetID, unrestrictedBooltargetId) ...
            uaWeightsAICMMA(estCoefs, y, covars, 'mma');
methodsArray{4}.shortName = 'mma';        
methodsArray{4}.longName = 'MMA';


methodsArray{5}.weightFunction = ...
    @(y, covars, estCoefs, estCovars, gradientEstimateTarget, ...
      targetID, unrestrictedBooltargetId) ...
        uaWeightsOptimal(estCoefs, estCovars, gradientEstimateTarget, ...
                         targetID, unrestrictedBooltargetId);
methodsArray{5}.shortName = 'opt';                   
methodsArray{5}.longName = 'Generic Optimal'; 
                    

%% Simulation
% Main simulation file
uaSimulateExp

%% Export plots
% The uaExportFigures script exports relative MSE plots 
% uaExportFigures