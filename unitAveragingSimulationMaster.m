clear variables
close all

%% Setup

plotQuietly = 0; % whether plots should be shown or saved without showing
% Define parameters used in all simulations here

% Data dimensions
Nvalues= [10, 25, 50]; % values of N to compare, must be arranged in ascending order
T = 60  ;



% For Dynamic model: select default mean and variance parameters for
% coefficients
lambdaMean = 0;
meanBeta = 1; % mean for beta coefficient
varianceBeta = 1; % variance for beta 
varNoiseVar = 1; % variance of heteroskedasticity
varianceX = 1; % to match the factor

% Number of replications
numReplications = 2500;

% Select thresholds for large N regime to kick in. We consider two values
% for comparison
i0_1 = 11;
i0_2 = 21;

localAsy = 1;

% Select a seed to drawing coefficients
seedCoefficients = 1; 


% Parameters
mu{1} = @(theta, x, y) theta(1, :);
mu{2} = @(theta, x, y) theta(1, :).*y(end, 1)+theta(2, :); % set x to 1
mu{3} = @(theta, x, y) theta(2, :)./(1-theta(1, :));
mu{4} = @(theta, x, y) theta(2, :); % beta
% Gradients
D{1} = @(theta, x, y) [1;0]; % lambda_1
D{2} = @(theta, x, y) [y(end, 1); 1]; % forecasts
D{3} = @(theta, x, y) [theta(2, 1)/(1-theta(1, 1))^2; 1/(1-theta(1, 1))] ; % long-run effect
D{4} = @(theta, x, y) [0; 1];
% Descriptions
des{1} = "\mu(\theta_1) = \lambda_1";
des{2} = "\mu(\theta_1) = E(y_{T+1}|y_T, x_T=1)";
des{3} = "\mu(\theta_1) = {\beta_1}/{1-\lambda_1}";
des{4} = "\mu(\theta_1) = \beta_1";
% Saving names, to generate figures
desName{1} = "lambda";
desName{2} = "forecast";
desName{3} = "longRun";
desName{4} = "beta";

%% Simulation
uaSimulate2

%% Export images
uaExportFigures

 
