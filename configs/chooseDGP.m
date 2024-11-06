% ===========================================================
% File: chooseDGP.m
% Description: This script controls the choice of DGPs for the simulations
% ===========================================================
%
% Project Name: Unit Averaging for Heterogeneous Panels
% Developed by: Christian Brownlees, Vladislav Morozov
 
% Implemented DGPs: "unimodal", "bimodal", "bimodal_close"
% If a list of DGPs is supplied, all of them will be ran in turn.

% Set DGPs for evaluating the mean squared error of unit averaging
simulationSettings = ["unimodal", "bimodal", "bimodal_close"];

% Set DGPs for evaluating the coverage of confidence intervals
simulationSettingsCI = ["ci_unimodal","ci_bimodal"];
