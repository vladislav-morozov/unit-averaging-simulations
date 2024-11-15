% ===========================================================
% Project Name: Unit Averaging for Heterogeneous Panels
% File: main.m
% Developed by: Christian Brownlees, Vladislav Morozov
% Contact: Vladislav Morozov (vladislav.morozov (at) barcelonagse.eu)
%
%
% Description:
% This script produces the simulation results reported in the paper and the
% Online Appendix.
%
%
% Link to paper: https://arxiv.org/abs/2210.14205
% Link to the appendix: 
% https://vladislav-morozov.github.io/files/1_unitAveragingSupplement.pdf
% 
%
% Acknowledgments:
% This project uses functions by the MATLAB community: 
% - ' ' by  
% 
%
% Last Modified: 2024-03-17
% Developed with MATLAB Version: R2022b
% ===========================================================

%% Initialization

% Clear the workspace
clc
clear variables
close all

% Load in folders with required functions and scripts
addpath('averagingFunctions') 
addpath(genpath('utilities')) % other utilities
addpath('configs') % Parameters of simulation
addpath('simulationScriptsFunctions') 

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

% Set data generating processes for the simulation
chooseDGP

% Set parameters of interest
chooseTargetParameters

% Load in blueprints of unit averaging approaches
% Note: optimal schemes implemented in createOptimalSchemes function
averagingApproachBlueprints
 
%% Simulation 1: performance of unit averaging approaches in terms of MSE

% Number of replications
numReplicationsMSE = 33;

for designID = 1:length(simulationSettingsMSE)
    % Extract current design name
    simulationSetting = simulationSettingsMSE(designID);
    
    % Set parameters
    setParameters 
    
    % Run simulation
    simulateMSE
    
    % Export plots 
    uaExportFigures
end 

%% Simulation 2: Performance of confidence intervas

% Number of replications
numReplicationsCI = 500;

for designID = 1:length(simulationSettingsCI)
    % Extract current design name
    simulationSetting = simulationSettingsCI(designID);
    
    % Set parameters
    setParameters
    simulateCI;
end 