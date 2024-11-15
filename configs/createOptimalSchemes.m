function optimalSchemes = ...
    createOptimalSchemes(thetaHat, thetaTrue, thetaLabels, ...
    averagingMode, averagingIncludeBool)
% createOptimalSchemes Implements the fixed-N and various large-N optimal 
% schemes
%
% Args:
%     thetaHat (matrix): A kxN matrix of coefficient estimates where 
%       columns index units.
%     thetaTrue (matrix): A kxN matrix of true coefficient vectors where 
%       columns index units.
%     thetaLabels (vector): An N-vector of integers where each component 
%       indicates the mixture component of the corresponding coefficient 
%       vector.
%     averagingMode (string): 'all' or 'firstOnly'. Determines if all units
%       are to serve as targets or only the first one.
%     averagingIncludeBool (vector): A boolean vector where True in the kth 
%       position means the kth approach is returned.
%
% Returns:
%     optimalSchemes (cell array): A cell array of structs where each
%       struct describes an averaging approach. Each approach is 
%       characterized by:
%         .shortName (string): Short name of the approach.
%         .longName (string): Long name of the approach.
%         .unrestrictedArray (cell array): An array of functions where each
%             function describes how to determine the unrestricted units 
%             for averaging the jth unit with a given scheme.
%          .(see averagingApproachBlueprints for further fields)
%
% Example:
%     thetaHat = rand(5, 10);  % 5 coefficients for 10 units
%     thetaTrue = rand(5, 10);
%     thetaLabels = randi([1, 3], 10, 1);
%     averagingMode = 'all';
%     averagingIncludeBool = true(1, 10);
%     optimalSchemes = optimalSchemes(thetaHat, thetaTrue, thetaLabels,...
%       averagingMode, averagingIncludeBool);
%

% Determine the number of targets
if averagingMode == "all"
    numTargets = size(thetaHat, 2);
else
    numTargets = 1;  % Only one unit will serve as the target
end

% Extract number of units
numUnits = size(thetaHat, 2);

% Initialize cell array to hold approaches
optimalSchemes = {};

%%%%%%%%%%%%%%%
%%% Fixed-N %%%
%%%%%%%%%%%%%%%

optimalSchemes{1}.shortName = 'unrestr';
optimalSchemes{1}.longName = 'Fixed-N';
for targetID = 1:numTargets
    optimalSchemes{1}.unrestrictedArray{targetID} = ...
        @(weightVector, paramFun) true(numUnits, 1);
end
optimalSchemes{1}.color =   [0, 0, 255]/255;  
optimalSchemes{1}.colorBW = [0.3, 0.3, 0.3];
optimalSchemes{1}.lineStyle = '-';
optimalSchemes{1}.marker = 'o';
optimalSchemes{1}.markerSize = 4;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Large-N: Focus oracles %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 10 closest units
optimalSchemes{2}.shortName = 'focus_oracle_10';
optimalSchemes{2}.longName = 'Large-N (10 closest focus)';
for targetID = 1:numTargets
    optimalSchemes{2}.unrestrictedArray{targetID} = ...
        @(weightVector, paramFun) ...
        oracleFocusSimilarUnits(thetaTrue, paramFun, targetID, 10/numUnits);
end
optimalSchemes{2}.color =   [15, 0, 124]/255;   
optimalSchemes{2}.colorBW = [0.6, 0.6, 0.6];
optimalSchemes{2}.lineStyle = ':'; 
optimalSchemes{2}.marker = '<';
optimalSchemes{2}.markerSize = 4;

% 25 closest units
optimalSchemes{3}.shortName = 'focus_oracle_25';
optimalSchemes{3}.longName = 'Large-N (25 closest focus)';
for targetID = 1:numTargets
    optimalSchemes{3}.unrestrictedArray{targetID} = ...
        @(weightVector, paramFun) ...
        oracleFocusSimilarUnits(thetaTrue, paramFun, targetID, 25/numUnits);
end
optimalSchemes{3}.color =   [150, 140, 100]/255;   
optimalSchemes{3}.colorBW = [0.6, 0.6, 0.6];
optimalSchemes{3}.lineStyle = ':'; 
optimalSchemes{3}.marker = 'v';
optimalSchemes{3}.markerSize = 4;

% 10% closest units
optimalSchemes{4}.shortName = 'focus_oracle_10_pct';
optimalSchemes{4}.longName = 'Large-N (10\% closest focus)';
for targetID = 1:numTargets
    optimalSchemes{4}.unrestrictedArray{targetID} = ...
        @(weightVector, paramFun) ...
        oracleFocusSimilarUnits(thetaTrue, paramFun, targetID, 0.1);
end
optimalSchemes{4}.color = [200, 200, 30]/255;   
optimalSchemes{4}.colorBW = [0.6, 0.6, 0.6];
optimalSchemes{4}.lineStyle = ':'; 
optimalSchemes{4}.marker = '*';
optimalSchemes{4}.markerSize = 4;

% 25% closest units
optimalSchemes{5}.shortName = 'focus_oracle_25_pct';
optimalSchemes{5}.longName = 'Large-N (25\% closest focus)';
for targetID = 1:numTargets
    optimalSchemes{5}.unrestrictedArray{targetID} = ...
        @(weightVector, paramFun) ...
        oracleFocusSimilarUnits(thetaTrue, paramFun, targetID, 0.25);
end
optimalSchemes{5}.color =   [1, 135, 232]/255; 
optimalSchemes{5}.colorBW = [0.6, 0.6, 0.6];
optimalSchemes{5}.lineStyle = ':'; 
optimalSchemes{5}.marker = 'o';
optimalSchemes{5}.markerSize = 4;

% 50% closest units
optimalSchemes{6}.shortName = 'focus_oracle_50_pct';
optimalSchemes{6}.longName = 'Large-N (50\% closest focus)';
for targetID = 1:numTargets
    optimalSchemes{6}.unrestrictedArray{targetID} = ...
        @(weightVector, paramFun) ...
        oracleFocusSimilarUnits(thetaTrue, paramFun, targetID, 0.5);
end
optimalSchemes{6}.color =   [0, 0, 0]/255;   
optimalSchemes{6}.colorBW = [0.6, 0.6, 0.6];
optimalSchemes{6}.lineStyle = ':'; 
optimalSchemes{6}.marker = 'square';
optimalSchemes{6}.markerSize = 4;


%%%%%%%%%%%%%%%%%%%%%%
%%% Large-N: Stein %%%
%%%%%%%%%%%%%%%%%%%%%%
 
optimalSchemes{7}.shortName = 'stein';
optimalSchemes{7}.longName = 'Large-N (Stein)';
for targetID = 1:numTargets
    optimalSchemes{7}.unrestrictedArray{targetID} = ...
        @(weightVector, paramFun)  ((1:numUnits)==targetID)';
end
optimalSchemes{7}.color =  [150, 140, 100]/255;   	
optimalSchemes{7}.colorBW = [0.5, 0.5, 0.5];
optimalSchemes{7}.lineStyle = '--'; 
optimalSchemes{7}.marker = 'pentagram';
optimalSchemes{7}.markerSize = 3;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Large-N: Focus cluster %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 

% 2 clusters
optimalSchemes{8}.shortName = 'focus_cluster_2';
optimalSchemes{8}.longName = 'Large-N (2 focus clusters)';
numClusters = 2;
for targetID = 1:numTargets 
    optimalSchemes{8}.unrestrictedArray{targetID} = ...
        @(weightVector, paramFun) ...
        focusClusterUnits(thetaTrue, paramFun, targetID, numClusters);
end
optimalSchemes{8}.color = [150, 140, 100]/255;   
optimalSchemes{8}.colorBW = [1, 1, 1];  
optimalSchemes{8}.lineStyle = ':'; 
optimalSchemes{8}.marker = '*';
optimalSchemes{8}.markerSize = 4;

% Main text: 4 clusters
optimalSchemes{9}.shortName = 'focus_cluster_4';
optimalSchemes{9}.longName = 'Large-N (4 focus clusters)';
numClusters = 4;
for targetID = 1:numTargets 
    optimalSchemes{9}.unrestrictedArray{targetID} = ...
        @(weightVector, paramFun) ...
        focusClusterUnits(thetaTrue, paramFun, targetID, numClusters);
end
optimalSchemes{9}.color =  [1, 135, 232]/255;  
optimalSchemes{9}.colorBW = [0.33, 0.33, 0.33];
optimalSchemes{9}.lineStyle = '--'; 
optimalSchemes{9}.marker = 'square';
optimalSchemes{9}.markerSize = 4;
optimalSchemes{9}.markerSize = 4;

% 8 clusters
optimalSchemes{10}.shortName = 'focus_cluster_8';
optimalSchemes{10}.longName = 'Large-N (8 focus clusters)';
numClusters = 8;
for targetID = 1:numTargets 
    optimalSchemes{10}.unrestrictedArray{targetID} = ...
        @(weightVector, paramFun) ...
        focusClusterUnits(thetaTrue, paramFun, targetID, numClusters);
end
optimalSchemes{10}.color =  [15, 0, 124]/255;  
optimalSchemes{10}.colorBW = [0.33, 0.33, 0.33];
optimalSchemes{10}.lineStyle = '--'; 
optimalSchemes{10}.marker = 'square';
optimalSchemes{10}.markerSize = 4;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Fixed-N: Top units %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Top weights of fixed-N criterion (fixed-N must also be available)

% Top 10 weights of fixed-N criterion  
optimalSchemes{11}.shortName = 'top_10';
optimalSchemes{11}.longName = 'Large-N (top 10 units)';
for targetID = 1:numTargets
    optimalSchemes{11}.unrestrictedArray{targetID} = ...
        @(weightVector, paramFun) boolTopCoords(weightVector, ...
        targetID, 10/numUnits);
end
optimalSchemes{11}.color =   [200, 200, 30]/255;   
optimalSchemes{11}.colorBW = [0, 0, 0];
optimalSchemes{11}.lineStyle = ':'; 
optimalSchemes{11}.marker = 'x';
optimalSchemes{11}.markerSize = 8;

% Top 25 weights of fixed-N criterion  
optimalSchemes{12}.shortName = 'top_25';
optimalSchemes{12}.longName = 'Large-N (top 25 units)';
for targetID = 1:numTargets
    optimalSchemes{12}.unrestrictedArray{targetID} = ...
        @(weightVector, paramFun) boolTopCoords(weightVector, ...
        targetID, 25/numUnits);
end
optimalSchemes{12}.color =   [150, 140, 100]/255;   
optimalSchemes{12}.lineStyle = ':'; 
optimalSchemes{12}.marker = '>';
optimalSchemes{12}.markerSize = 8;

% Top 10%
optimalSchemes{13}.shortName = 'top_10_pct';
optimalSchemes{13}.longName = 'Large-N (top 10\% units)';
for targetID = 1:numTargets
    optimalSchemes{13}.unrestrictedArray{targetID} = ...
        @(weightVector, paramFun) boolTopCoords(weightVector, ...
        targetID, 0.1);
end
optimalSchemes{13}.color =   [170, 120, 10]/250;  
optimalSchemes{13}.colorBW = [0, 0, 0];
optimalSchemes{13}.lineStyle = '--'; 
optimalSchemes{13}.marker = '<';
optimalSchemes{13}.markerSize = 8;

% Top 25%
optimalSchemes{14}.shortName = 'top_25_pct';
optimalSchemes{14}.longName = 'Large-N (top 25\% units)';
for targetID = 1:numTargets
    optimalSchemes{14}.unrestrictedArray{targetID} = ...
        @(weightVector, paramFun) boolTopCoords(weightVector, ...
        targetID, 0.25);
end
optimalSchemes{14}.color =   [1, 135, 232]/255; 
optimalSchemes{14}.colorBW = [0, 0, 0];
optimalSchemes{14}.lineStyle = '--'; 
optimalSchemes{14}.marker = 'o';
optimalSchemes{14}.markerSize = 8;

% Top 50%
optimalSchemes{15}.shortName = 'top_50_pct';
optimalSchemes{15}.longName = 'Large-N (top 50\% units)';
for targetID = 1:numTargets
    optimalSchemes{15}.unrestrictedArray{targetID} = ...
        @(weightVector, paramFun) boolTopCoords(weightVector, ...
        targetID, 0.5);
end
optimalSchemes{15}.color =   [15, 0, 124]/255;   
optimalSchemes{15}.colorBW = [0, 0, 0];
optimalSchemes{15}.lineStyle = '--'; 
optimalSchemes{15}.marker = 'square';
optimalSchemes{15}.markerSize = 8;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Large-N: Random restrictions %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Random 10 units
optimalSchemes{16}.shortName = 'random_10';
optimalSchemes{16}.longName = 'Large-N (random 10)';
for targetID = 1:numTargets
    randomBool = false(numUnits, 1);
    randomBool(datasample(1:numUnits, min(numUnits, 10), ...
        'Replace',false)) = true;
    randomBool(targetID) = true;
    optimalSchemes{16}.unrestrictedArray{targetID} = ...
        @(weightVector, paramFun) randomBool;
end
optimalSchemes{16}.color =   [150, 140, 100]/255;   
optimalSchemes{16}.lineStyle = ':'; 
optimalSchemes{16}.marker = '>';
optimalSchemes{16}.markerSize = 4;

% Random 20 units
optimalSchemes{17}.shortName = 'random_20';
optimalSchemes{17}.longName = 'Large-N (random 20)';
for targetID = 1:numTargets
    randomBool = false(numUnits, 1);
    randomBool(datasample(1:numUnits, min(numUnits, 20), ...
        'Replace',false)) ...
        = true;
    randomBool(targetID) = true;
    optimalSchemes{17}.unrestrictedArray{targetID} = ...
        @(weightVector, paramFun) randomBool;
end
optimalSchemes{17}.color =   [1, 135, 232]/255;  
optimalSchemes{17}.lineStyle = ':'; 
optimalSchemes{17}.marker = '<';
optimalSchemes{17}.markerSize = 4; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Large-N: Oracle cluster %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

optimalSchemes{18}.shortName = 'oracleClasses';
optimalSchemes{18}.longName = 'Large-N (true class labels)';
for targetID = 1:numTargets
    currentCluster = thetaLabels(targetID);
    optimalSchemes{18}.unrestrictedArray{targetID} = ...
        @(weightVector, paramFun) thetaLabels == currentCluster;
end
optimalSchemes{18}.color =   [38, 185, 159]/255;   
optimalSchemes{18}.lineStyle = ':'; 
optimalSchemes{18}.colorBW =   [0.55, 0.55, 0.55];
optimalSchemes{18}.marker = '+';
optimalSchemes{18}.markerSize = 7;



% Anti-oracle labels
optimalSchemes{19}.shortName = 'antiOracleClasses';
optimalSchemes{19}.longName = 'Large-N (opposite class labels)';
for targetID = 1:numTargets
    currentCluster = thetaLabels(targetID);
    optimalSchemes{19}.unrestrictedArray{targetID} = ...
        @(weightVector, paramFun) thetaLabels ~= currentCluster;
end
optimalSchemes{19}.color =   [38, 100, 200]/255;   
optimalSchemes{19}.lineStyle = ':'; 
optimalSchemes{19}.marker = '^';

% Extract only the desired averaging schemes
optimalSchemes = optimalSchemes(averagingIncludeBool);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Large-N: Coefficient Oracles %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% 10 closest
optimalSchemes{20}.shortName = 'oracle_similar_10';
optimalSchemes{20}.longName = 'Large-N (10 closest coef)';
for targetID = 1:numTargets
    dists = sum((thetaTrue-thetaTrue(:, 1)).^2);
    [~, I] = sort(dists); 
    mostSimilarIdx = I(1:10);
    similarBool = false(numUnits, 1);
    similarBool(mostSimilarIdx) = true;
    optimalSchemes{20}.unrestrictedArray{targetID} = ...
        @(weightVector, paramFun) similarBool;
end
optimalSchemes{20}.color =   [15, 0, 124]/255;  %  dark blue
optimalSchemes{20}.colorBW = [0.6, 0.6, 0.6];
optimalSchemes{20}.lineStyle = ':'; 
optimalSchemes{20}.marker = '<';
optimalSchemes{20}.markerSize = 4;

% 25 closest
optimalSchemes{21}.shortName = 'oracle_similar_25';
optimalSchemes{21}.longName = 'Large-N (25 closest coef)';
for targetID = 1:numTargets
    dists = sum((thetaTrue-thetaTrue(:, 1)).^2);
    [~, I] = sort(dists); 
    mostSimilarIdx = I(1:25);
    similarBool = false(numUnits, 1);
    similarBool(mostSimilarIdx) = true;
    optimalSchemes{21}.unrestrictedArray{targetID} = ...
        @(weightVector, paramFun) similarBool;
end
optimalSchemes{21}.color =   [150, 140, 100]/255;   
optimalSchemes{21}.colorBW = [0.6, 0.6, 0.6];
optimalSchemes{21}.lineStyle = ':'; 
optimalSchemes{21}.marker = 'v';
optimalSchemes{21}.markerSize = 4;

% 10% closest
optimalSchemes{22}.shortName = 'oracle_similar_10_pct';
optimalSchemes{22}.longName = 'Large-N (10\% closest coef)';
for targetID = 1:numTargets
    dists = sum((thetaTrue-thetaTrue(:, 1)).^2);
    [~, I] = sort(dists); 
    mostSimilarIdx = I(1:ceil(0.1*numUnits));
    similarBool = false(numUnits, 1);
    similarBool(mostSimilarIdx) = true;
    optimalSchemes{22}.unrestrictedArray{targetID} = ...
        @(weightVector, paramFun) similarBool;
end
optimalSchemes{22}.color = [200, 200, 30]/255;   
optimalSchemes{22}.colorBW = [0.6, 0.6, 0.6];
optimalSchemes{22}.lineStyle = ':'; 
optimalSchemes{22}.marker = '*';
optimalSchemes{22}.markerSize = 4;

% 25% closest
optimalSchemes{23}.shortName = 'oracle_similar_25_pct';
optimalSchemes{23}.longName = 'Large-N (25\% closest coef)';
for targetID = 1:numTargets
    dists = sum((thetaTrue-thetaTrue(:, 1)).^2);
    [~, I] = sort(dists); 
    mostSimilarIdx = I(1:ceil(0.25*numUnits));
    similarBool = false(numUnits, 1);
    similarBool(mostSimilarIdx) = true;
    optimalSchemes{23}.unrestrictedArray{targetID} = ...
        @(weightVector, paramFun) similarBool;
end
optimalSchemes{23}.color =   [1, 135, 232]/255;  
optimalSchemes{23}.colorBW = [0.6, 0.6, 0.6];
optimalSchemes{23}.lineStyle = ':'; 
optimalSchemes{23}.marker = 'o';
optimalSchemes{23}.markerSize = 4;

% 50% closest
optimalSchemes{24}.shortName = 'oracle_similar_50_pct';
optimalSchemes{24}.longName = 'Large-N (50\% closest coef)';
for targetID = 1:numTargets
    dists = sum((thetaTrue-thetaTrue(:, 1)).^2);
    [~, I] = sort(dists); 
    mostSimilarIdx = I(1:ceil(0.5*numUnits));
    similarBool = false(numUnits, 1);
    similarBool(mostSimilarIdx) = true;
    optimalSchemes{24}.unrestrictedArray{targetID} = ...
        @(weightVector, paramFun) similarBool;
end
optimalSchemes{24}.color =   [0, 0, 0]/255;   
optimalSchemes{24}.colorBW = [0.6, 0.6, 0.6];
optimalSchemes{24}.lineStyle = ':'; 
optimalSchemes{24}.marker = 'square';
optimalSchemes{24}.markerSize = 4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Large-N: Coef clusters %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 2 clusters
coefClusters = kmeans(thetaHat', 2);  
optimalSchemes{25}.shortName = 'cluster_coef_2';
optimalSchemes{25}.longName = 'Large-N (2 coef clusters)';
for targetID = 1:numTargets
    currentCluster = coefClusters(targetID);
    optimalSchemes{25}.unrestrictedArray{targetID} = ...
        @(weightVector, paramFun) coefClusters == currentCluster;
end
optimalSchemes{25}.color = [150, 140, 100]/255;   
optimalSchemes{25}.colorBW = [1, 1, 1];  
optimalSchemes{25}.lineStyle = ':'; 
optimalSchemes{25}.marker = '*';
optimalSchemes{25}.markerSize = 4;

% Main text: 4 clusters
coefClusters = kmeans(thetaHat', 4); 
optimalSchemes{26}.shortName = 'cluster_coef_4';
optimalSchemes{26}.longName = 'Large-N (4 coef clusters)';
for targetID = 1:numTargets
    currentCluster = coefClusters(targetID);
    optimalSchemes{26}.unrestrictedArray{targetID} = ...
        @(weightVector, paramFun) coefClusters == currentCluster;
end
optimalSchemes{26}.color =  [1, 135, 232]/255;  
optimalSchemes{26}.colorBW = [0.33, 0.33, 0.33];
optimalSchemes{26}.lineStyle = '--'; 
optimalSchemes{26}.marker = 'square';
optimalSchemes{26}.markerSize = 4;
optimalSchemes{26}.markerSize = 4;

% 8 clusters
coefClusters = kmeans(thetaHat', 8);  
optimalSchemes{27}.shortName = 'cluster_coef_8';
optimalSchemes{27}.longName = 'Large-N (8 coef clusters)';
for targetID = 1:numTargets
    currentCluster = coefClusters(targetID);
    optimalSchemes{27}.unrestrictedArray{targetID} = ...
        @(weightVector, paramFun) coefClusters == currentCluster;
end
optimalSchemes{27}.color =  [15, 0, 124]/255;   
optimalSchemes{27}.colorBW = [0.33, 0.33, 0.33];
optimalSchemes{27}.lineStyle = '--'; 
optimalSchemes{27}.marker = 'square';
optimalSchemes{27}.markerSize = 4;
end


function targetBool = boolTopCoords(weightVector, targetIdx, topShare)
% boolTopCoords Returns a Boolean vector indicating the top share of
% weights.
%
% Args:
%     weightVector (vector): N-vector of averaging weights.
%     targetIdx (int): Position of the target vector.
%     topShare (double): Between 0 and 1, indicating the share of units to 
%       return.
%
% Returns:
%     targetBool (logical vector): N-vector where the top 
%         topShare*100 percent of weights are marked as true, including the
%         target position.

% Extract the dimension of the problem
numUnits = length(weightVector);
% Determine the number of top values to extract
numberTop = ceil(topShare * numUnits);
% Identify the indices of the top values
[~, topIds] = maxk(weightVector, numberTop);
% Create output Boolean vector
targetBool = false(numUnits, 1);
targetBool(topIds) = true;
targetBool(targetIdx) = true;

end


function similarBool = ...
    oracleFocusSimilarUnits(thetaTrue, paramFun, targetID, share)
% oracleFocusSimilarUnits Identifies units closest to the target unit.
%
% Args:
%     thetaTrue (matrix): kxN sample of true values.
%     paramFun (function handle): Function to apply to thetaTrue.
%     targetID (int): Position of the target unit.
%     share (double): Between 0 and 1, indicating the share of units to 
%         flag based on closeness.
%
% Returns:
%     similarBool (logical vector): N-vector where the top share*100  
%         percent of units closest to the target unit are marked as true.
%

% Extract size
[~, numUnits] = size(thetaTrue);

% Apply the function to the sample
paramSample = paramFun(thetaTrue);

% Evaluate distances to the target unit
dists = abs(paramSample - paramSample(targetID));

% Sort distances and identify the closest units
[~, I] = sort(dists);
mostSimilarIdx = I(1:ceil(share * numUnits));

% Create output Boolean vector
similarBool = false(numUnits, 1);
similarBool(mostSimilarIdx) = true;

end


function clusterBool = ...
    focusClusterUnits(thetaVector, paramFun, targetID, numClusters)
% focusClusterUnits Clusters units and returns units in the same cluster as
% the target unit.
%
% Args:
%     thetaVector (matrix): kxN sample of true values.
%     paramFun (function handle): Function to apply to thetaVector.
%     targetID (int): Position of the target unit.
%     numClusters (int): Number of clusters to use.
%
% Returns:
%     clusterBool (logical vector): N-vector where the units in the same  
%         cluster as the target unit are marked as true.

% Apply the function to the sample
paramSample = paramFun(thetaVector);

% Perform clustering
focusClusters = kmeans(paramSample', numClusters);

% Identify the cluster of the target unit
currentCluster = focusClusters(targetID);

% Create output Boolean vector
clusterBool = focusClusters == currentCluster;
end
