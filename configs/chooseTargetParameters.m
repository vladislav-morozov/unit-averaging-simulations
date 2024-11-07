% ===========================================================
% File: chooseTargetParameters.m
% Description: This script defines target parameters for simulation
%              estimation in the paper.
% ===========================================================
%
% Project Name: Unit Averaging for Heterogeneous Panels
% Authors: Christian Brownlees, Vladislav Morozov
%
% Model: Heterogeneous panel ARX(1) model
%   y_{it} = theta_{i1} y_{it-1} + theta_{i2} x_{it} + u_{it}
%
% Parameters of Interest:
%   This script defines target parameters (mu(theta)) for estimation:
%   1. theta_{i1} - Coefficient on lagged variable y_{it-1}
%   2. theta_{i2} - Coefficient on exogenous covariate x_{it}
%   3. One-step ahead forecast based on most recent y and predefined 
%      xForecast
%
% Parameter Array:
%   - Each target parameter is defined as an element in `paramArray`.
%   - Each parameter is a struct with fields:
%       - mu: The target function mu(theta)
%       - gradient: Gradient of mu(theta) with respect to coefficients
%       - plotDescr: Description of parameter for plot labels
%       - saveName: Identifier used in filenames for saved plots
%
% Note: Parameters are stored in `paramArray`, with each parameter accessed
%       as `paramArray{i}`. The script iterates over the array length to 
%       access all defined parameters.
%
% ===========================================================
%% Parameter 1: Coefficient on Lag (theta_{i1})
% mu(theta) function: Returns theta_{i1}, the coefficient on the lag
paramArray{1}.mu = @(theta, x, y) theta(1, :);

% Gradient of mu with respect to theta
paramArray{1}.gradient = @(theta, x, y) [1; 0];

% Plotting descriptor for theta_{i1}
paramArray{1}.plotDescr = "$\mu(\theta_1) = \lambda_1$";

% Save name identifier for this parameter
paramArray{1}.saveName = "lambda";

%% Parameter 2: Coefficient on Exogenous Variable (theta_{i2})
% mu(theta) function: Returns theta_{i2}, the coefficient on x_{it}
paramArray{2}.mu = @(theta, x, y) theta(2, :);

% Gradient of mu with respect to theta
paramArray{2}.gradient = @(theta, x, y) [0; 1];

% Plotting descriptor for theta_{i2}
paramArray{2}.plotDescr = "$\mu(\theta_2) = \beta_1$";

% Save name identifier for this parameter
paramArray{2}.saveName = "beta";

%% Parameter 3: One-Step Ahead Forecast
% Value of exogenous variable for forecasting
xForecast = 1;

% mu(theta) function: Forecast based on the last observed y and xForecast
paramArray{3}.mu = @(theta, x, y) theta(1, :) .* y(end, 1) + ...
                                  theta(2, :) * xForecast;

% Gradient of forecast function with respect to theta
paramArray{3}.gradient = @(theta, x, y) [y(end, 1); 1];

% Plotting descriptor for one-step ahead forecast
paramArray{3}.plotDescr = "$\mu(\theta) = E(y_{1,T+1}|y_T, x_{1,T+1}=1)$";

% Save name identifier for forecast parameter
paramArray{3}.saveName = "forecast";
