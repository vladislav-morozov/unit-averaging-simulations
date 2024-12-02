% ===========================================================
% File: setApproachesToPlot.m
% Description: This script creates bundles of averaging approaches that
%              will be plotted together
% ===========================================================
%
% Project Name: Unit Averaging for Heterogeneous Panels
% Developed by: Christian Brownlees, Vladislav Morozov
%
% This script determines which collections of unit averaging approaches get
% plotted. A separate plot is generated for each collection.
%
% Each collection is defined based on three variables:
%   - Vector of strings containing the names of the approaches to plots.
%   - A struct indicating whether any of the approaches should be renamed
%     in the legend relative to their default names
%   - A string that will form part of the plot name.
% ===========================================================

%% Create bundles of approaches
% Specified using short names of averaging approaches

% Plots reported in the main text
approachesToPlotMSEBW = ["ind", "mg", "aic", ...
    "unrestr", "oracle_similar_10_pct","stein", ...
    "top_10_pct"];
% Previous main plot: excluded for being too crowded
approachesToPlotMSEPaper = ["ind", "mg", "aic", ...
    "unrestr", "oracle_similar_10_pct", "stein", ...
    "top_10_pct","cluster_coef_4"];
% Oracle unrestricted units
approachesToPlotFocusOracleCome = ...
    ["focus_oracle_10", "focus_oracle_25", ...
    "focus_oracle_10_pct","focus_oracle_25_pct", "focus_oracle_50_pct"];
% Empirical similarity of parameters
approachesToPlotSimilComp = ...
    ["oracle_similar_10", "oracle_similar_25" , ...
    "oracle_similar_10_pct", "oracle_similar_25_pct", ...
    "oracle_similar_50_pct"];
% Clustering based on coefficients
approachesToPlotClusterComp = ...
    ["cluster_coef_2", "cluster_coef_4", "cluster_coef_8"];
% Clustering based on the parameters
approachesToPlotFocusClusterComp = ...
    ["focus_cluster_2", "focus_cluster_4", "focus_cluster_8"];

% If the approach is multimodal, add information on true classes
if coefApproach ~= "unimodal"
    approachesToPlotMSEPaper = ...
        [approachesToPlotMSEPaper, "oracleClasses"];
    approachesToPlotSimilComp = ...
        [approachesToPlotSimilComp, "oracleClasses"];
    approachesToPlotClusterComp = ...
        [approachesToPlotClusterComp, "oracleClasses"];
end

% Design-independent approaches
approachesToPlotTopComp = ...
    ["top_10", "top_25", "top_10_pct", "top_25_pct", "top_50_pct"];
approachesToPlotRandomComp = ...
    ["random_10", "random_20"];
approachesToPlotFirstWeight = ...
    ["unrestr", "top", "oracleSimilarity", "stein"];
approachesAlt = [ "unrestr","focus_oracle_10_pct",...
    "top_10_pct", "focus_cluster_4"];
approachesAlt2 =  [ "unrestr","focus_oracle_50_pct",...
    "top_10_pct", "focus_cluster_2"];
approachesToPlotAnimated = ["ind", "mg", "unrestr"];

%% Patching names of plotted approaches

% --- Main text renaming ---
% Simpler names for large-N approaches
approachPaperRenames.oracle_similar_10_pct = "Large-N (most similar)";
approachPaperRenames.oracle_similar_10 = "Large-N (most similar)";
approachPaperRenames.top_10_pct = "Large-N (top units)";
approachPaperRenames.top_10 = "Large-N (top units)";
approachPaperRenames.cluster_coef_4 = "Large-N (cluster coefs)";

% --- Online Appendix renaming ---
% Add a * to the names of the approaches reported in the main text
approachOARenamesTop.top_10_pct = ...
    "Large-N (top 10 units, *)";
approachOARenamesSimil.oracle_similar_10_pct = ...
    "Large-N (10 most similar, *)";
approachOARenamesCluster.cluster_coef_4 = ...
    "Large-N (4 coef clusters)";

% Empty renames: no renames needed
approachOARenamesRandom = struct();
approachOARenamesFocusCluster = struct();
approachOARenamesOracleFocus = struct();

% Renames for alternative possibilities of picking large-N approaches
approachGridRenames.focus_oracle_10_pct = "Large-N (most similar)";
approachGridRenames.focus_oracle_25_pct = "Large-N (most similar)";
approachGridRenames.focus_oracle_50_pct = "Large-N (most similar)";
approachGridRenames.focus_cluster_2 = "Large-N (cluster coefs)";
approachGridRenames.focus_cluster_4 = "Large-N (cluster coefs)";
approachGridRenames.focus_cluster_8 = "Large-N (cluster coefs)";
approachGridRenames.top_10_pct = "Large-N (top units)";

%% Combine line bundles together with renames, attach a saving name

lineSets{1} = approachesToPlotMSEPaper;
nameSets{1} = approachPaperRenames;
plotName{1} = "color";

lineSets{2} = approachesToPlotTopComp;
nameSets{2} = approachOARenamesTop;
plotName{2}= "comp_top";

lineSets{3} = approachesToPlotSimilComp;
nameSets{3} = approachOARenamesSimil;
plotName{3}= "comp_simil";

lineSets{4} = approachesToPlotClusterComp;
nameSets{4} = approachOARenamesCluster;
plotName{4}= "comp_cluster";

lineSets{5} = approachesToPlotRandomComp;
nameSets{5} = approachOARenamesRandom;
plotName{5}= "comp_random";

lineSets{6} = approachesToPlotFocusClusterComp;
nameSets{6} = approachOARenamesFocusCluster;
plotName{6}= "comp_focus_cluster";

lineSets{7} = approachesToPlotFocusOracleCome;
nameSets{7} = approachOARenamesOracleFocus;
plotName{7} = "comp_focus_oracle";

lineSets{8} = approachesAlt2;
nameSets{8} = approachGridRenames;
plotName{8} = "diff_lines";


%% Extract the approaches that were computed with the current settings

optimalSchemes = ...
    createOptimalSchemes(randn(2, N), randn(2, N), randn(N, 1), ...
    'firstOnly', averagingIncludeBool);
allMethodsArray = ...
    addOptimalToMethodsStruct(methodsArray, optimalSchemes);