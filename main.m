% ===========================================================
% Project Name: Unit Averaging for Heterogeneous Panels
% File: main.m
% Developed by: Christian Brownlees, Vladislav Morozov
% Contact: Vladislav Morozov (vladislav.morozov (at) barcelonagse.eu)
% ===========================================================
%
% Description: 
%   This script reproduces the simulation results presented in the 
%   associated research paper and its Online Appendix. It evaluates the 
%   performance of various unit averaging approaches in terms of:
%       - Mean Squared Error (MSE)
%       - Bias
%       - Variance
%       - Weights assigned to units
%       - Probabilites of units being unrestricted in large-N approaches
%       - Coverage and length properties of confidence intervals.
%
% Usage: 
%   To reproduce the full set of simulation results, execute this script 
%   directly. The script will generate all necessary outputs, including 
%   figures and tables, which are automatically saved to the 'results' 
%   folder.
%
% Output: 
%   - Figures: Saved in PDF and PNG formats;
%   - Tables: Exported as LaTeX files.
%
% For simulation design and overall background on unit averaging:
%   - Link to paper: 
%       arxiv.org/abs/2210.14205
%   - Link to Online Appendix: 
%       vladislav-morozov.github.io/assets/files/1_unitAveragingSupplement.pdf
%
% Software Requirements:
%   - MATLAB (tested on R2022b and R2024b). 
%   - Toolboxes: 
%       - Parallel Computing Toolbox (used for simulation speed-up). 
%   - Additional File Exchange Dependencies:
%       - `table2latex` (for LaTeX table generation). 
%       - `tight_subplot` (for subplot layout optimization). 
%         Both dependencies are provided in the replication package.
%
% Notes:
%   - Ensure all dependencies are in the MATLAB path before execution.
%   - The simulation is computationally intensive and may benefit from 
%     running on a machine with multiple CPU cores.
% 
% Reference runtimes:
%   - Around 28 hours on a Intel i9-14900K CPU.
%
% Developed with MATLAB Versions: R2022b to R2024b
% ===========================================================

%% Initialization

% Clear the workspace
clc
clear variables
close all

% Load in folders with required functions and scripts 
addpath(genpath('src'))         % configs and functions
addpath(genpath('scripts'))     % simulation and export scripts

%% Set the parameters of the simulation

% Whether plots should be shown or saved without showing
plotQuietly = false; 

% Whether the weights and the IDs of the unrestricted units should be saved
% for each sample.  
saveWeights = true;
saveUnrestricted = false; 
 
% Select a seed to draw data
seedCoefficients = 1; 

% Set number of samples to draw
numReplicationsMSE = 10000;      % For MSE simulation
numReplicationsCI = 500;        % For confidence interval simulations

% Set data generating processes for the simulation
chooseDGP

% Set parameters of interest
chooseTargetParameters

% Load in blueprints of unit averaging approaches
% Note: optimal schemes implemented in createOptimalSchemes function
averagingApproachBlueprints
 
%% Simulation 1: performance of unit averaging approaches in terms of MSE
% Evaluates MSE, bias, variance, weights and unrestricted units 

for designID = 1:length(simulationSettingsMSE)
    % Extract current design name
    simulationSetting = simulationSettingsMSE(designID);
    
    % Set parameters
    setParameters 
    
    % Run simulation
    simulateMSE
    
    % Export plots 
    exportFigures
end 

%% Simulation 2: Performance of confidence intervas
% Evaluates coverage and length properties of confidence intervals 

for designID = 1:length(simulationSettingsCI)
    % Extract current design name
    simulationSetting = simulationSettingsCI(designID);
    
    % Set parameters
    setParameters
    
    % Run confidence interval simulations
    simulateCI;

    % Export tables
    exportCoverageLengthTables
end 