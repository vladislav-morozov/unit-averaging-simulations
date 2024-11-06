% ===========================================================
% File: averagingApproachBlueprints.m
% Description: This script implements blueprints for averaging approaches
% ===========================================================
%
% Project Name: Unit Averaging for Heterogeneous Panels
% Developed by: Christian Brownlees, Vladislav Morozovn.
%
%  
% IMPLEMENTATION NOTES:
% - The generic optimal scheme is used as a blueprint for other optimal
%   schemes. As some large-N schemes depend on the data, these must be
%   created for each given sample. The different schemes are implemented by
%   the createOptimalSchemes function
% - Each method is defined by a struct. The fields are the weight function,
%   names to use for plotting, and plotting parameters (colors, line
%   styles)
%   Implemented before improvement in performance of Matlab classes

%% Individual estimator
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


%%  Mean group
methodsArray{2}.weightFunction = ...
    @(y, covars, estCoefs, estCovars, gradientEstimateTarget, ...
      targetID, unrestrictedBooltargetId) ...
            ones(size(estCoefs, 2), 1)/size(estCoefs, 2);
methodsArray{2}.shortName = 'mg';        
methodsArray{2}.longName = 'Mean Group';
methodsArray{2}.color = [0.7, 0.7, 0.7]; 
methodsArray{2}.colorBW = [0.79, 0.79, 0.75];
methodsArray{2}.lineStyle = '-.';
methodsArray{2}.marker = 'none';
methodsArray{2}.markerSize = 6;


%% Exponential AIC weights
methodsArray{3}.weightFunction = ...
    @(y, covars, estCoefs, estCovars, gradientEstimateTarget, ...
      targetID, unrestrictedBooltargetId) ...
            uaWeightsAICMMA(estCoefs, y, covars, 'aic');
methodsArray{3}.shortName = 'aic';        
methodsArray{3}.longName = 'AIC';
methodsArray{3}.color =  [0.44, 0.44, 0.4];  
methodsArray{3}.colorBW = [0.5, 0.5, 0.5];
methodsArray{3}.lineStyle = '-.';
methodsArray{3}.marker = '+';
methodsArray{3}.markerSize = 8;


%% Mallows averaging
methodsArray{4}.weightFunction = ...
    @(y, covars, estCoefs, estCovars, gradientEstimateTarget, ...
      targetID, unrestrictedBooltargetId) ...
            uaWeightsAICMMA(estCoefs, y, covars, 'mma');
methodsArray{4}.shortName = 'mma';        
methodsArray{4}.longName = 'MMA';
methodsArray{4}.color =   [ 55, 130, 125]/255; 
methodsArray{4}.colorBW = [0.01, 0.01, 0.01]; % unplotted in BW
methodsArray{4}.lineStyle = '-.';
methodsArray{4}.marker = '.';
methodsArray{4}.markerSize = 6;


%% Generic optimal approach
methodsArray{5}.weightFunction = ...
    @(y, covars, estCoefs, estCovars, gradientEstimateTarget, ...
      targetID, unrestrictedBooltargetId) ...
        uaWeightsOptimal(estCoefs, estCovars, gradientEstimateTarget, ...
                         targetID, unrestrictedBooltargetId);
methodsArray{5}.shortName = 'opt';                   
methodsArray{5}.longName = 'Generic Optimal'; 
% The colors and styles should be set for individual units, see
methodsArray{5}.color =   [0, 0, 0];
methodsArray{5}.colorBW = [0, 0, 0];
methodsArray{5}.lineStyle = '--';
methodsArray{5}.marker = 'x';
methodsArray{5}.markerSize = 6;