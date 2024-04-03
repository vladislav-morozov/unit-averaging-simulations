%% (Mostly) common parameters

% variance of heteroskedasticity
varNoiseVar = 1;
% Variance of covariates
varianceX = 1;

% Default set: different for local and CI simulations
averagingIncludeBool = logical([1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1]);

%% Setting-Specific Parameters
if simulationSetting == "unimodal"
    coefApproach = 'unimodal';
 
    % Parameter range
    theta1Range = 0.2:.04:0.8;
    
    % Data dimensions
    valuesN= [50, 150, 450];  
    valuesT = [30, 60, 180, 600] ;
   
    % Plotting parameters
    yMinRelMSEPlot = 0.65;
    yMaxRelMSEPlot = 1.51;
elseif simulationSetting == "bimodal"
    coefApproach = 'bimodal';
 
    % Parameter range
    theta1Range = 0.1:0.05:0.8;
    
    % Data dimensions
    valuesN= [50, 150, 450]; 
    valuesT = [30, 60, 180, 600] ;

    % Plotting parameters
    yMinRelMSEPlot = 0.72;
    yMaxRelMSEPlot = 1.51;
elseif simulationSetting == "local"
    % Approach of v1 of the paper -- simulations in the local framework
    
    % Coefficient DGP
    coefApproach = "local";
    % variance of heteroskedasticity
    varNoiseVar = 1;
    varianceX = 1;
    
    % Parameter range
    theta1Range = 0:0.04:0.99;
    
    % Data dimensions
    valuesN= [50, 150, 250]; % values of N to compare, must be arranged in ascending order
    valuesT = 60 ;
    
    yMinRelMSEPlot = 0.72;
    yMaxRelMSEPlot = 1.51;
    yMinAbsMSEPlot = 0;
    yMaxAbsMSEPlot = 0.015;
    
    % Methods to plot: match the methods reported in the original paper
    averagingIncludeBool = logical([1, 0, 0, 0, 1, 1, 1]);
    
    
elseif simulationSetting == "ci_unimodal"
    
    
    % Unimodal setting for CI evaluation
    coefApproach = 'unimodal';
 
    
    valuesN= [50, 150, 450]; % values of N to compare, must be arranged in ascending order
    valuesT = [60, 180] ;
    
    
    % Parameter range
    theta1Range = linspace(0.2, 0.8, 5);
    
    numBootstrapSamples= 500;
    
    % Only include fixed-N weights (unrestr)
    averagingIncludeBool = logical([1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]);
    
elseif simulationSetting == "ci_bimodal"
    
    
    % Unimodal setting for CI evaluation
    coefApproach = 'bimodal';
    
    
    valuesN= [50, 150, 450]; % values of N to compare, must be arranged in ascending order
    valuesT = [60, 180] ;
    
    
    % Parameter range
    theta1Range = linspace(0.2, 0.8, 5);
    
    numBootstrapSamples= 500;
    
    % Only include fixed-N weights (unrestr)
    averagingIncludeBool = logical([1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]);
end




%% Plotting parameters
% Set thegrid for evaluating weight behavior
thetaGrid = theta1Range;

approachesToPlotMSE = ["ind", "mg", "aic", "mma", "unrestr",...
    "cluster_coef", "top", "stein", "random10pct", "oracleSimilarity" ];
approachesToPlotFirstWeight = ["unrestr", "top", "oracleSimilarity", "stein"];
%% Create design-specific folder if it does not exist
warning('off', 'MATLAB:MKDIR:DirectoryExists');
outputFolderName = 'Outputs/'+simulationSetting;
figureFolderName  = 'Figures/'+simulationSetting;
mkdir(figureFolderName);
mkdir(outputFolderName);