% ===========================================================
% File: chooseTargetParameters.m
% Description: This script controls the choice of target parameters for the
% simulations of the paper
% ===========================================================
%
% Project Name: Unit Averaging for Heterogeneous Panels
% Developed by: Christian Brownlees, Vladislav Morozov
% MATLAB Version: R2022b
%
%
% The model is a heterogeneous panel ARX(1) model of the form:
%   y_{it} = theta_{i1} y_{it-1} + theta_{i2} x_{it} + u_{it}
%
%
% PARAMETERS CONSIDERED:
% The methodology applies to smooth transformation mu(theta)
% This script implements 3 parameters of interest (mu):
% 1. theta_{i1} -- coefficient on the lag
% 2. theta_{i2} -- coefficient on the exogenous variable 
% 3. one-step ahead out-of-sample forecast. Depends on the last value of y
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

 

% Coefficient on lag
paramArray{1}.mu = @(theta, x, y) theta(1, :); 
paramArray{1}.gradient = @(theta, x, y) [1;0]; 
paramArray{1}.plotDescr = "$\mu(\theta_1) = \lambda_1$";
paramArray{1}.saveName = "lambda";

% Coefficient on exogenous variable
paramArray{2}.mu = @(theta, x, y) theta(2, :);
paramArray{2}.gradient = @(theta, x, y) [0; 1]; 
paramArray{2}.plotDescr = "$\mu(\theta_1) = \beta_1$";
paramArray{2}.saveName = "beta";

% Forecast
xForecast = 1;
paramArray{3}.mu =  @(theta, x, y) ...
    theta(1, :).*y(end, 1)+theta(2, :)*xForecast; 
paramArray{3}.gradient = @(theta, x, y) [y(end, 1); 1];  
paramArray{3}.plotDescr = "$\mu(\theta_1) = E(y_{1T+1}|y_T, x_{1T+1}=1)$";
paramArray{3}.saveName = "forecast";
