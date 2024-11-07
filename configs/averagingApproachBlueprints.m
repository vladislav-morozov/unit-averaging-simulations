% ===========================================================
% File: averagingApproachBlueprints.m
% Description: This script defines various unit averaging approaches
%              for heterogeneous panel data.
% ===========================================================
%
% Project Name: Unit Averaging for Heterogeneous Panels
% Developed by: Christian Brownlees, Vladislav Morozov
%
% Implementation Notes:
%   - This script sets up a series of averaging methods using the 
%     generic optimal scheme as a basis for customization.
%   - Some large-N schemes require re-creation for each sample due to 
%     data dependencies and are implemented by the `createOptimalSchemes` 
%     function.
%   - Each method is represented by a struct with fields:
%       - weightFunction: Function for computing weights
%       - shortName: Short label for plots
%       - longName: Full name for legend/plotting
%       - color, colorBW: Colors for color and black-and-white plots
%       - lineStyle, marker, markerSize: Plotting style parameters
%
% ===========================================================

%% Averaging Method 1: Individual Estimator
% Weight function for individual estimators (selects individual unit)
methodsArray{1}.weightFunction = ...
    @(y, covars, estCoefs, estCovars, gradientEstimateTarget, ...
      targetID, unrestrictedBooltargetId) ...
            double(1:size(estCoefs, 2) == targetID)';

% Plotting names for the individual estimator
methodsArray{1}.shortName = 'ind';
methodsArray{1}.longName = 'Individual';

% Plotting parameters for individual estimator
methodsArray{1}.color = [0.01, 0.01, 0.05];
methodsArray{1}.colorBW = [0.01, 0.01, 0.01];
methodsArray{1}.lineStyle = '-';
methodsArray{1}.marker = 'none';
methodsArray{1}.markerSize = 6;

%% Averaging Method 2: Mean Group
% Weight function for mean group estimator (uniform weights)
methodsArray{2}.weightFunction = ...
    @(y, covars, estCoefs, estCovars, gradientEstimateTarget, ...
      targetID, unrestrictedBooltargetId) ...
            ones(size(estCoefs, 2), 1) / size(estCoefs, 2);

% Plotting names for the mean group estimator
methodsArray{2}.shortName = 'mg';
methodsArray{2}.longName = 'Mean Group';

% Plotting parameters for mean group estimator
methodsArray{2}.color = [0.7, 0.7, 0.7];
methodsArray{2}.colorBW = [0.79, 0.79, 0.75];
methodsArray{2}.lineStyle = '-.';
methodsArray{2}.marker = 'none';
methodsArray{2}.markerSize = 6;

%% Averaging Method 3: Exponential AIC Weights
% Weight function based on AIC (Akaike Information Criterion) weights
methodsArray{3}.weightFunction = ...
    @(y, covars, estCoefs, estCovars, gradientEstimateTarget, ...
      targetID, unrestrictedBooltargetId) ...
            uaWeightsAICMMA(estCoefs, y, covars, 'aic');

% Plotting names for AIC-based method
methodsArray{3}.shortName = 'aic';
methodsArray{3}.longName = 'AIC';

% Plotting parameters for AIC weights
methodsArray{3}.color = [0.44, 0.44, 0.4];
methodsArray{3}.colorBW = [0.5, 0.5, 0.5];
methodsArray{3}.lineStyle = '-.';
methodsArray{3}.marker = '+';
methodsArray{3}.markerSize = 8;

%% Averaging Method 4: Mallows Averaging (MMA)
% Weight function based on Mallows Model Averaging (MMA)
methodsArray{4}.weightFunction = ...
    @(y, covars, estCoefs, estCovars, gradientEstimateTarget, ...
      targetID, unrestrictedBooltargetId) ...
            uaWeightsAICMMA(estCoefs, y, covars, 'mma');

% Plotting names for Mallows averaging
methodsArray{4}.shortName = 'mma';
methodsArray{4}.longName = 'MMA';

% Plotting parameters for Mallows averaging
methodsArray{4}.color = [55, 130, 125] / 255;
methodsArray{4}.colorBW = [0.01, 0.01, 0.01]; % Not plotted in BW
methodsArray{4}.lineStyle = '-.';
methodsArray{4}.marker = '.';
methodsArray{4}.markerSize = 6;

%% Averaging Method 5: Generic Optimal Approach
% Weight function for generic optimal averaging approach
methodsArray{5}.weightFunction = ...
    @(y, covars, estCoefs, estCovars, gradientEstimateTarget, ...
      targetID, unrestrictedBooltargetId) ...
        uaWeightsOptimal(estCoefs, estCovars, gradientEstimateTarget, ...
                         targetID, unrestrictedBooltargetId);

% Plotting names for generic optimal approach
methodsArray{5}.shortName = 'opt';
methodsArray{5}.longName = 'Generic Optimal';

% Plotting parameters for generic optimal approach
methodsArray{5}.color = [0, 0, 0];
methodsArray{5}.colorBW = [0, 0, 0];
methodsArray{5}.lineStyle = '--';
methodsArray{5}.marker = 'x';
methodsArray{5}.markerSize = 6;
