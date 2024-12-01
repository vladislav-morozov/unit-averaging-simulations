% ===========================================================
% File: chooseGidPlots.m
% Description: This script creates defines the exported grid plots  
% ===========================================================
%
% Project Name: Unit Averaging for Heterogeneous Panels
% Developed by: Christian Brownlees, Vladislav Morozov
%
% This script describes the various 2-dimensional grip plots. Plot
% descriptions are saved in the cell array gridPlots.
%
% Each plot is described by a struct with the following fields
% - data (cell array): 
%       cell array of tables indexed by (N, T). Each element is a table 
%       with columns indexing unit averaging approaches.
% - approachesToPlot (vector): 
%       vector of strings containing the names of the approaches to plot. 
%       This attribute is modified at plot creation time with the current 
%       approach bundle. For creation, use the dummy list. 
% - logScale (Boolean): 
%       whether a log scale should be used for the data.
% - cmapN1, cmapN2 (int): 
%       number of colors to use in the colormaps
% - caxisLims (vector): 
%       limits for the color axis
% - suptitle (string): 
%       Suptitle for the plot.
% - plotSavingName (string): 
%       Used as part of the file name when saving the plot.
% ===========================================================

%% Average weight as a function of own and target values

gridPlots{1}.data = weightsTablesNT;
gridPlots{1}.approachesToPlot = approachesAlt;
gridPlots{1}.logScale = true;
gridPlots{1}.cmapN1 = 150;
gridPlots{1}.cmapN2 = 45;
gridPlots{1}.cmapFlip = true;
gridPlots{1}.caxisLims = [0.001, 1];
gridPlots{1}.plotSavingName = "weights_2d";
gridPlots{1}.suptitle = ...
    "Average weight of a unit with $\lambda_{alt}$ in estimating $\lambda_1$";

%% Probability of being unrestricted as a function of own and target values

gridPlots{2}.data = unitsUnrestrNT; 
gridPlots{2}.approachesToPlot = approachesAlt; 
gridPlots{2}.logScale = false;
gridPlots{2}.cmapN1 = 150;
gridPlots{2}.cmapN2 = 20;
gridPlots{2}.cmapFlip = true;
gridPlots{2}.caxisLims = [0, 1];
gridPlots{2}.plotSavingName = "unrestr_2d";
gridPlots{2}.suptitle = ...
    "Probability of a unit with $\lambda_{alt}$ being unrestricted in estimating $\lambda_1$";