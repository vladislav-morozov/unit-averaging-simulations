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
addpath('utilities')
addpath('configs') 

%% Set the parameters of the simulation
% Whether plots should be shown or saved without showing
plotQuietly = false; 

% Whether the weights and the IDs of the unrestricted units should be saved
% for each sample.  
saveWeights = true;
saveUnrestricted = false; 
 

% Select a seed to drawing coefficients
seedCoefficients = 1; 



% Set data generating processes for the simulation
chooseDGP

% Set parameters of interest
chooseTargetParameters

% Create unit averaging schemes

 
%% Averaging schemes
%
% NOTE: the generic optimal scheme is used as a blueprint for other optimal
% schemes. As some large-N schemes depend on the data, these must be created
% for each given sample. The different schemes are implemented by the
% uaAddOptimalToMethodsStruct function

methodsArray{1}.weightFunction = ...
    @(y, covars, estCoefs, estCovars, gradientEstimateTarget, ...
      targetID, unrestrictedBooltargetId) ...
            double(1:size(estCoefs, 2) == targetID)';
methodsArray{1}.shortName = 'ind';        
methodsArray{1}.longName = 'Individual';
methodsArray{1}.color = [0.01, 0.01, 0.05]; 
methodsArray{1}.colorBW = [0.01, 0.01, 0.01];  
methodsArray{1}.lineStyle = '-';
methodsArray{1}.marker = 'none';
methodsArray{1}.markerSize = 6;



methodsArray{2}.weightFunction = ...
    @(y, covars, estCoefs, estCovars, gradientEstimateTarget, ...
      targetID, unrestrictedBooltargetId) ...
            ones(size(estCoefs, 2), 1)/size(estCoefs, 2);
methodsArray{2}.shortName = 'mg';        
methodsArray{2}.longName = 'Mean Group';
methodsArray{2}.color = [0.7, 0.7, 0.7]; %[141, 111, 100]/255; % brown
methodsArray{2}.colorBW = [0.79, 0.79, 0.75];
methodsArray{2}.lineStyle = '-.';
methodsArray{2}.marker = 'none';
methodsArray{2}.markerSize = 6;

methodsArray{3}.weightFunction = ...
    @(y, covars, estCoefs, estCovars, gradientEstimateTarget, ...
      targetID, unrestrictedBooltargetId) ...
            uaWeightsAICMMA(estCoefs, y, covars, 'aic');
methodsArray{3}.shortName = 'aic';        
methodsArray{3}.longName = 'AIC';
methodsArray{3}.color =  [0.44, 0.44, 0.4]; % [ 106, 130, 68]/255;
methodsArray{3}.colorBW = [0.5, 0.5, 0.5];
methodsArray{3}.lineStyle = '-.';
methodsArray{3}.marker = '+';
methodsArray{3}.markerSize = 8;

methodsArray{4}.weightFunction = ...
    @(y, covars, estCoefs, estCovars, gradientEstimateTarget, ...
      targetID, unrestrictedBooltargetId) ...
            uaWeightsAICMMA(estCoefs, y, covars, 'mma');
methodsArray{4}.shortName = 'mma';        
methodsArray{4}.longName = 'MMA';
methodsArray{4}.color =   [ 55, 130, 125]/255; % greenish-teal
methodsArray{4}.colorBW = [0.01, 0.01, 0.01]; % unplotted in BW
methodsArray{4}.lineStyle = '-.';
methodsArray{4}.marker = '.';
methodsArray{4}.markerSize = 6;


methodsArray{5}.weightFunction = ...
    @(y, covars, estCoefs, estCovars, gradientEstimateTarget, ...
      targetID, unrestrictedBooltargetId) ...
        uaWeightsOptimal(estCoefs, estCovars, gradientEstimateTarget, ...
                         targetID, unrestrictedBooltargetId);
methodsArray{5}.shortName = 'opt';                   
methodsArray{5}.longName = 'Generic Optimal'; 
% The colors and styles should be set for individual units
methodsArray{5}.color =   [0, 0, 0];
methodsArray{5}.colorBW = [0, 0, 0];
methodsArray{5}.lineStyle = '--';
methodsArray{5}.marker = 'x';
methodsArray{5}.markerSize = 6;

%% MSE Simulation 

% Number of replications
numReplications = 3333;

for designID = 1:length(simulationSettings)
    % Extract current design name
    simulationSetting = simulationSettings(designID);
    
    % Close all figures
    close all
    
    % Set parameters
    uaSetParameters 
    
    % Run simulation
    uaSimulate
    
    % Export plots 
    uaExportFigures
end 

%% Coverage simulation

% Number of replications
numReplications = 500;

for designID = 1:length(simulationSettingsCI)
    % Extract current design name
    simulationSetting = simulationSettingsCI(designID);
    
    % Set parameters
    uaSetParameters
    uaSimulateCI;
end 