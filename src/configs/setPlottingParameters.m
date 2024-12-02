% ===========================================================
% File: setPlottingParameters.m
% Description: This script sets parameters related to figure creation
% ===========================================================
%
% Project Name: Unit Averaging for Heterogeneous Panels
% Developed by: Christian Brownlees, Vladislav Morozov
% ===========================================================

% Set text interpreters to latex
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter', 'latex');

% Scalar parameters
plotLineThickness= 1.5;     % Line thickness
plotWLines = 1000;          % Width of figure in pixels
plotHLines = 800;           % Height in figure pixels

% Grid of points for line plots. Markers will be placed at positions
% theta1Range only
gridMultiplier = 3;  
thetaGridMSE = ...
    linspace(min(theta1Range), max(theta1Range), ...
             gridMultiplier*length(theta1Range));

% Axis plotting limits
yLimsRelMSE.lambda = [0.5, 1.3];     % for the AR(1) parameter
yLimsRelMSE.beta = [0.9, 1.2];       % for the coefficient on x_{it}
yLimsRelMSE.forecast = [0.8, 1.3];   % for the forecast parameters