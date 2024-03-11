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
numReplications = 200;


%% DGP Parameters
% Data dimensions
valuesN= [10, 25, 50, 100, 150, 250]; % values of N to compare, must be arranged in ascending order
valuesT = 60 ;

% For Dynamic model: select default mean and variance parameters for
% coefficients
lambdaMean = 0;
meanBeta = 1; % mean for beta coefficient
varianceBeta = 1; % variance for beta 
varNoiseVar = 1; % variance of heteroskedasticity
varianceX = 1; % to match the factor

%% Large-N specifications
% Implementation notes


% Select thresholds for large N regime to kick in. We consider two values
% for comparison
i0_1 = 16;
i0_2 = 51;

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
% 2. one-step ahead out-of-sample forecast. Depends on the last value of y
%    and a value of xForecast specified in the definition
% 3. long-run effect of a one-unit shift in x
% 4. theta_{i2} -- coefficient on the exogenous variable
%
%
% IMPLEMENTATION NOTES: 
% -- the function mu takes both the coefficients and the
%    full data vectors y and x. 
% -- the gradient of mu must be supplied explicitly as a separation
%    function D
% -- the plotting functions will use the description and name of the
%    parameters to automatically generate labels and file names of the
%    plots
% -- the script detects the length of the mu array and loops through
%    all the parameters

% Parameters
mu{1} = @(theta, x, y) theta(1, :); 
xForecast = 1;
mu{2} = @(theta, x, y) theta(1, :).*y(end, 1)+theta(2, :)*xForecast; 
mu{3} = @(theta, x, y) theta(2, :)./(1-theta(1, :));
mu{4} = @(theta, x, y) theta(2, :); 
% Gradients
D{1} = @(theta, x, y) [1;0]; 
D{2} = @(theta, x, y) [y(end, 1); 1]; 
D{3} = @(theta, x, y) [theta(2, 1)/(1-theta(1, 1))^2; 1/(1-theta(1, 1))];
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
% Main simulation file
uaSimulate

%% Export plots
% The uaExportFigures script exports relative MSE plots 
uaExportFigures