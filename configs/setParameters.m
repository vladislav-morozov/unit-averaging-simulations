% ===========================================================
% File: setParameters.m
% Description: This script sets parameters for various simulation designs
% ===========================================================
%
% Project Name: Unit Averaging for Heterogeneous Panels
% Developed by: Christian Brownlees, Vladislav Morozov
%
% Model Description: Heterogeneous panel ARX(1) model
%   y_{it} = theta_{i1} y_{it-1} + theta_{i2} x_{it} + u_{it}
%
% Implemented DGPs for coefficient:
%   "unimodal", "bimodal", "bimodal_close", "ci_unimodal", "ci_bimodal",
%   "local"
%
% Common Parameters:
%   - varNoiseVar: Variance of error term u_{it}
%   - varianceX: Variance of covariate x_{it}
%   - averagingIncludeBool: Boolean array to control inclusion of unit 
%       averaging approaches
%
% Setting-Specific Parameters:
%   Parameters vary by simulation setting (e.g., theta1Range, sample sizes)
%
% ===========================================================
%% Common Parameters
% Variance of the noise term u_{it}
varNoiseVar = 1;

% Variance of exogenous covariates x_{it}
varianceX = 1;

% Boolean vector: determines which unit averaging approaches to include
% By default, all implemented approaches are considered, and modified
% in specific settings (e.g., CI evaluation for confidence intervals).
averagingIncludeBool = true(1, 19);

%% Setting-Specific Parameters
% Configure parameters based on the selected simulation setting
switch simulationSetting
    case "unimodal"
        % Unimodal setting
        coefApproach = 'unimodal';
        theta1Range = 0.2:0.04:0.8;     % Range for theta_{i1}
        valuesN = [50, 150, 450];       % Sample sizes
        valuesT = [30, 60, 180, 600];   % Time points

    case "bimodal"
        % Bimodal setting
        coefApproach = 'bimodal';
        theta1Range = 0.2:0.05:0.8;
        valuesN = [50, 150, 450];
        valuesT = [30, 60, 180, 600];

    case "bimodal_close"
        % Bimodal setting with closer modes
        coefApproach = 'bimodal_close';
        theta1Range = 0.2:0.05:0.8;
        valuesN = [50, 150, 450];
        valuesT = 60;

    case "local"
        % Local framework simulation (as in v1 of the paper)
        coefApproach = "local";
        theta1Range = 0:0.04:0.99;
        valuesN = [50, 150, 250];
        valuesT = 60;

        % Match methods to those reported in the original paper
        averagingIncludeBool = logical([1, 0, 0, 0, 1, 1, 1]);

    case "ci_unimodal"
        % Unimodal setting for confidence interval (CI) evaluation
        coefApproach = 'unimodal';
        valuesN = [50, 150, 450];
        valuesT = [60, 180];
        theta1Range = linspace(0.2, 0.8, 5);
        numBootstrapSamples = 500; % number of samples to draw for CIs

        % Include only fixed-N weights
        averagingIncludeBool = logical([1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]);

    case "ci_bimodal"
        % Bimodal setting for CI evaluation
        coefApproach = 'bimodal';
        valuesN = [50, 150, 450];
        valuesT = [60, 180];
        theta1Range = linspace(0.2, 0.8, 5);
        numBootstrapSamples = 500;

        % Include only fixed-N weights
        averagingIncludeBool = logical([1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]);

    otherwise
        error('Unknown simulation setting: %s', simulationSetting);
end

%% Plotting Parameters
% Grid for plotting results (based on theta1Range from the selected
% setting)
thetaGrid = theta1Range;

%% Output Directory Setup
% Create directories for saving simulation results and figures if they
% don't exist
outputFolderName = fullfile('Outputs', simulationSetting);
figureFolderName = fullfile('Figures', simulationSetting);

% Suppress directory creation warnings if directories already exist
warning('off', 'MATLAB:MKDIR:DirectoryExists');
mkdir(figureFolderName);
mkdir(outputFolderName);
warning('on', 'MATLAB:MKDIR:DirectoryExists'); % Re-enable warnings
