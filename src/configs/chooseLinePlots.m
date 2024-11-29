% ===========================================================
% File: chooseLinePlots.m
% Description: This script creates defines the line plots exported
% ===========================================================
%
% Project Name: Unit Averaging for Heterogeneous Panels
% Developed by: Christian Brownlees, Vladislav Morozov
%
% This script describes the various 1-dimensional line plots. Plot
% descriptions are saved in the cell array linePlots.
%
% Each plot is described by a struct with the following fields
% - data (cell array): 
%       cell array of tables indexed by (N, T). Each element is a table 
%       with columns indexing unit averaging approaches.
% - approachesToPlot (vector): 
%       vector of strings containing the names of the approaches to plot. 
%       This attribute is modified at plot creation time with the current 
%       approach bundle. For creation, use the dummy list.
% - firstColumnPlot (integer): 
%       columns to the left of this one will be ignored. Used for defining 
%       plots in relation to the individual estimator (always in column 1).
% - dataTransform (function): 
%       which aggregation should be applied to the data field to produce 
%       current line.
% - yLim (vector): 
%       limits for the y-axis
% - relative (Boolean): 
%       if true, a horizontal line y=1 is drawn.
% - suptitle (string): 
%       Suptitle for the plot.
% - plotSavingName (string): 
%       Used as part of the file name when saving the plot.
% ===========================================================

% Dummy placeholder list of approaches
approachesToPlotDummy = "unrestr";

%% Mean squared error relative to the individual estimator

linePlots{1}.data = mseTablesNT;                         % data 
linePlots{1}.approachesToPlot = approachesToPlotDummy;   % approach names
linePlots{1}.firstColumnPlot = 2;                        % first column
linePlots{1}.dataTransform = ...                         % data transform
    @(dataTable, colIdx) ...
    dataTable{:, colIdx}./dataTable{:, 1};
linePlots{1}.yLims = @(dataTable, paramName) ...         % y-axis limits
    [yLimsRelMSE.(paramName)];
linePlots{1}.relative = true;                            % add y=1?
linePlots{1}.suptitle = @(plotDescr) ...                 % suptitle
    {'Averaging estimators, ' + ...
    plotDescr,  'Ratio of MSE to individual estimator'};
linePlots{1}.plotSavingName = 'relative_mse';            % saving name

%% Absolute MSE

linePlots{2}.data = mseTablesNT;
linePlots{2}.approachesToPlot = approachesToPlotDummy;
linePlots{2}.firstColumnPlot = 1;
linePlots{2}.dataTransform = ...
    @(dataTable, colIdx) dataTable{:, colIdx};
linePlots{2}.yLims = @(dataTable, paramName) ...
    [0, 1.5*max(dataTable.unrestr)];
linePlots{2}.relative = false;
linePlots{2}.suptitle = @(plotDescr) 'Averaging estimators, ' + ...
    plotDescr + ...
    ', MSE';
linePlots{2}.plotSavingName = 'absolute_mse';

%% Bias

linePlots{3}.data = biasTablesNT;
linePlots{3}.approachesToPlot = approachesToPlotDummy;
linePlots{3}.firstColumnPlot = 1;
linePlots{3}.dataTransform = ...
    @(dataTable, colIdx) (dataTable{:, colIdx});
linePlots{3}.yLims = @(dataTable, paramName) ...
    [1.5*min(dataTable.unrestr), 1.5*max(dataTable.unrestr)];
linePlots{3}.relative = false;
linePlots{3}.suptitle = @(plotDescr) 'Averaging estimators, ' + ...
    plotDescr + ...
    ', bias';
linePlots{3}.plotSavingName = 'bias';

%% Variance relative to the individual estimator

linePlots{4}.data = varTablesNT;
linePlots{4}.approachesToPlot = approachesToPlotDummy;
linePlots{4}.firstColumnPlot = 2;
linePlots{4}.dataTransform = ...
    @(dataTable, colIdx) (dataTable{:, colIdx})./(dataTable{:, 1});
linePlots{4}.yLims = @(dataTable, paramName) [0, ...
    1.5*max(  (dataTable.unrestr./(dataTable{:, 1})))];
linePlots{4}.relative = true;
linePlots{4}.suptitle = @(plotDescr) 'Averaging estimators, ' + ...
    plotDescr + ...
    ', relative variance';
linePlots{4}.plotSavingName = 'relative_variance';

%% Average weight of first unit

linePlots{5}.data = firstWeightNT;
linePlots{5}.approachesToPlot = approachesToPlotDummy;
linePlots{5}.firstColumnPlot = 1;
linePlots{5}.dataTransform = ...
    @(dataTable, colIdx) dataTable{:, colIdx};
linePlots{5}.yLims = @(dataTable, paramName) [0, 1];
linePlots{5}.relative = false;
linePlots{5}.suptitle = @(plotDescr) {'Minimum MSE estimators, ' + ...
    plotDescr, ...
    'Average weight of target unit'};
linePlots{5}.plotSavingName = 'average_first_weight';

%% Max difference in weights relative to fixed-N estimator 

linePlots{6}.data = maxDiffNT;
linePlots{6}.approachesToPlot = approachesToPlotDummy;
linePlots{6}.firstColumnPlot = 2;
linePlots{6}.dataTransform = ...
    @(dataTable, colIdx) dataTable{:, colIdx};
linePlots{6}.yLims = @(dataTable, paramName) [0, 0.52];
linePlots{6}.relative = false;
linePlots{6}.suptitle = @(plotDescr)...
    {'Average maximum difference in weights of unrestricted units',  ...
    'vs. fixed-N estimator, ' + ...
    plotDescr};
linePlots{6}.plotSavingName = 'max_diff';

%% Mass difference relative to fixed-N estimator 

linePlots{7}.data = massDiffNT;
linePlots{7}.approachesToPlot = approachesToPlotDummy;
linePlots{7}.firstColumnPlot = 2;
linePlots{7}.dataTransform = ...
    @(dataTable, colIdx) dataTable{:, colIdx};
linePlots{7}.yLims = @(dataTable, paramName) [-0.7, 0.2];
linePlots{7}.relative = false;
linePlots{7}.suptitle = @(plotDescr) ...
    {'Average difference in total mass of unrestricted units', ...
    'vs. fixed-N estimator, ' + ...
    plotDescr};
linePlots{7}.plotSavingName = 'mass_diff';