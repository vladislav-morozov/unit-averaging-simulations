
 

if simulationSetting == "unimodal"
    coefApproach = 'unimodal';
    % variance of heteroskedasticity
    varNoiseVar = 2; 
    varianceX = 2;  

    % Parameter range
    theta1Range = 0.2:.04:0.8;
    
    % Data dimensions
    valuesN= [50, 150, 450]; % values of N to compare, must be arranged in ascending order
    valuesT = [30, 60, 180, 600] ;

    averagingIncludeBool = logical([1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0]);
    yMinRelMSEPlot = 0.65;
    yMaxRelMSEPlot = 1.51;
elseif simulationSetting == "bimodal"
    coefApproach = 'bimodal';
    % variance of heteroskedasticity
    varNoiseVar = 2; 
    varianceX = 2;  

    % Parameter range
    theta1Range = 0.1:0.05:0.8;
    
    % Data dimensions
    valuesN= [50, 150, 450]; % values of N to compare, must be arranged in ascending order
    valuesT = [30, 60, 180, 600] ;

    % Schemes to include
    averagingIncludeBool = logical([1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0]);
    
    % Plotting parameters
    yMinRelMSEPlot = 0.72;
    yMaxRelMSEPlot = 1.51;
elseif simulationSetting == "local"
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
    
    % Methods to plot 
    averagingIncludeBool = logical([1, 0, 0, 0, 1, 1]);
end
%% Common parameters

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