function optimalSchemes = ...
    uaOptimalSchemes(thetaHat, thetaTrue, thetaLabels,...
    averagingMode, averagingIncludeBool)
% uaOptimalSchemes Implements the fixed-N and various large-N optimal
% schemes. 6 approaches are implemented currently
%
% Inputs: 1. thetaHat -- kxN matrix of coefficient estimates. Columns index
%            units
%         2. averagingMode -- 'all' or 'firstOnly'. Determines if all units
%            are to serve as targets or only the first one 
%         3.
%         4. 
%         5. averagingIncludeBool -- boolean 6-vector. True in kth position
%            means kth approach is returned

%
% Returns: optimalSchemes -- cell array of structs. Each struct describes a
%   an averaging approaches. An approach is characterized by 3 fields:
%   .shortName, .longName, and .unrestrictedArray. The latter is an array 
%   of functions. The jth function describes the way to determine the
%   unrestricted units for averaging the jth unit with a given scheme.


if averagingMode == "all"
    numTargets = size(thetaHat, 2);
else
    numTargets = 1;
end

numUnits =  size(thetaHat, 2);

% Do any sample-specific large-N computations
optimalSchemes = {};

%%%%%%%%%%%%%%%
%%% Fixed-N %%%
%%%%%%%%%%%%%%%

optimalSchemes{1}.shortName = 'unrestr';
optimalSchemes{1}.longName = 'Fixed-N';
for targetID = 1:numTargets
    optimalSchemes{1}.unrestrictedArray{targetID} = ...
        @(weightVector) true(numUnits, 1);
end
optimalSchemes{1}.color =   [ 132, 1, 3]/255;  % dark red
optimalSchemes{1}.colorBW = [0.3, 0.3, 0.3];
optimalSchemes{1}.lineStyle = '-';
optimalSchemes{1}.marker = 'o';
optimalSchemes{1}.markerSize = 4;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Similarity Oracles %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

% Main text one: 10 closest
optimalSchemes{2}.shortName = 'oracle_similar_10';
optimalSchemes{2}.longName = 'Large-N (10 most similar)';
for targetID = 1:numTargets
    dists = sum((thetaTrue-thetaTrue(:, 1)).^2);
    [~, I] = sort(dists); 
    mostSimilarIdx = I(1:10);
    similarBool = false(numUnits, 1);
    similarBool(mostSimilarIdx) = true;
    optimalSchemes{2}.unrestrictedArray{targetID} = ...
        @(weightVector) similarBool;
end
optimalSchemes{2}.color =   [0, 0, 255]/255;  %  most blue
optimalSchemes{2}.colorBW = [0.6, 0.6, 0.6];
optimalSchemes{2}.lineStyle = ':'; 
optimalSchemes{2}.marker = '<';
optimalSchemes{2}.markerSize = 4;

% 25 closest
optimalSchemes{3}.shortName = 'oracle_similar_25';
optimalSchemes{3}.longName = 'Large-N (25 most similar)';
for targetID = 1:numTargets
    dists = sum((thetaTrue-thetaTrue(:, 1)).^2);
    [~, I] = sort(dists); 
    mostSimilarIdx = I(1:10);
    similarBool = false(numUnits, 1);
    similarBool(mostSimilarIdx) = true;
    optimalSchemes{3}.unrestrictedArray{targetID} = ...
        @(weightVector) similarBool;
end
optimalSchemes{3}.color =   [0, 0, 255]/255;  %  most blue
optimalSchemes{3}.colorBW = [0.6, 0.6, 0.6];
optimalSchemes{3}.lineStyle = ':'; 
optimalSchemes{3}.marker = '<';
optimalSchemes{3}.markerSize = 4;

% 10% closest
optimalSchemes{4}.shortName = 'oracle_similar_10_pct';
optimalSchemes{4}.longName = 'Large-N (10\% most similar)';
for targetID = 1:numTargets
    dists = sum((thetaTrue-thetaTrue(:, 1)).^2);
    [~, I] = sort(dists); 
    mostSimilarIdx = I(1:ceil(0.1*numUnits));
    similarBool = false(numUnits, 1);
    similarBool(mostSimilarIdx) = true;
    optimalSchemes{4}.unrestrictedArray{targetID} = ...
        @(weightVector) similarBool;
end
optimalSchemes{4}.color =   [0, 0, 255]/255;  %  most blue
optimalSchemes{4}.colorBW = [0.6, 0.6, 0.6];
optimalSchemes{4}.lineStyle = ':'; 
optimalSchemes{4}.marker = '<';
optimalSchemes{4}.markerSize = 4;

% 25% closest
optimalSchemes{5}.shortName = 'oracle_similar_25_pct';
optimalSchemes{5}.longName = 'Large-N (25\% most similar)';
for targetID = 1:numTargets
    dists = sum((thetaTrue-thetaTrue(:, 1)).^2);
    [~, I] = sort(dists); 
    mostSimilarIdx = I(1:ceil(0.25*numUnits));
    similarBool = false(numUnits, 1);
    similarBool(mostSimilarIdx) = true;
    optimalSchemes{5}.unrestrictedArray{targetID} = ...
        @(weightVector) similarBool;
end
optimalSchemes{5}.color =   [0, 0, 255]/255;  %  most blue
optimalSchemes{5}.colorBW = [0.6, 0.6, 0.6];
optimalSchemes{5}.lineStyle = ':'; 
optimalSchemes{5}.marker = '<';
optimalSchemes{5}.markerSize = 4;

% 50% closest
optimalSchemes{6}.shortName = 'oracle_similar_50_pct';
optimalSchemes{6}.longName = 'Large-N (50\% most similar)';
for targetID = 1:numTargets
    dists = sum((thetaTrue-thetaTrue(:, 1)).^2);
    [~, I] = sort(dists); 
    mostSimilarIdx = I(1:ceil(0.5*numUnits));
    similarBool = false(numUnits, 1);
    similarBool(mostSimilarIdx) = true;
    optimalSchemes{6}.unrestrictedArray{targetID} = ...
        @(weightVector) similarBool;
end
optimalSchemes{6}.color =   [0, 0, 255]/255;  %  most blue
optimalSchemes{6}.colorBW = [0.6, 0.6, 0.6];
optimalSchemes{6}.lineStyle = ':'; 
optimalSchemes{6}.marker = '<';
optimalSchemes{6}.markerSize = 4;


%%%%%%%%%%%%%
%%% Stein %%%
%%%%%%%%%%%%%
 
optimalSchemes{7}.shortName = 'stein';
optimalSchemes{7}.longName = 'Large-N (Stein)';
for targetID = 1:numTargets
    optimalSchemes{7}.unrestrictedArray{targetID} = ...
        @(weightVector)  ((1:numUnits)==targetID)';
end
optimalSchemes{7}.color =   [179, 112, 79]/255;  % weird brown
optimalSchemes{7}.colorBW = [0.5, 0.5, 0.5];
optimalSchemes{7}.lineStyle = '--'; 
optimalSchemes{7}.marker = 'pentagram';
optimalSchemes{7}.markerSize = 3;


%%%%%%%%%%%%%%%%%%
%%% Clustering %%%
%%%%%%%%%%%%%%%%%%

% 2 clusters
coefClusters = kmeans(thetaHat', 2); % cluster coefs
optimalSchemes{8}.shortName = 'cluster_coef_2';
optimalSchemes{8}.longName = 'Large-N (2 coef clusters)';
for targetID = 1:numTargets
    currentCluster = coefClusters(targetID);
    optimalSchemes{8}.unrestrictedArray{targetID} = ...
        @(weightVector) coefClusters == currentCluster;
end
optimalSchemes{8}.color =  [90, 30, 236]/255;  % purple
optimalSchemes{8}.colorBW = [1, 1, 1]; % unplotted in BW
optimalSchemes{8}.lineStyle = ':'; 
optimalSchemes{8}.marker = '*';

% Main text: 4 clusters
coefClusters = kmeans(thetaHat', 4); % cluster coefs
optimalSchemes{9}.shortName = 'cluster_coef_4';
optimalSchemes{9}.longName = 'Large-N (4 coef clusters)';
for targetID = 1:numTargets
    currentCluster = coefClusters(targetID);
    optimalSchemes{9}.unrestrictedArray{targetID} = ...
        @(weightVector) coefClusters == currentCluster;
end
optimalSchemes{9}.color =  [218, 1, 136]/255;  % intense pink
optimalSchemes{9}.colorBW = [0.33, 0.33, 0.33];
optimalSchemes{9}.lineStyle = '--'; 
optimalSchemes{9}.marker = 'square';
optimalSchemes{9}.markerSize = 4;

% 8 clusters
coefClusters = kmeans(thetaHat', 8); % cluster coefs
optimalSchemes{10}.shortName = 'cluster_coef_8';
optimalSchemes{10}.longName = 'Large-N (8 coef clusters)';
for targetID = 1:numTargets
    currentCluster = coefClusters(targetID);
    optimalSchemes{10}.unrestrictedArray{targetID} = ...
        @(weightVector) coefClusters == currentCluster;
end
optimalSchemes{10}.color =  [218, 1, 136]/255;  % intense pink
optimalSchemes{10}.colorBW = [0.33, 0.33, 0.33];
optimalSchemes{10}.lineStyle = '--'; 
optimalSchemes{10}.marker = 'square';
optimalSchemes{10}.markerSize = 4;

%%%%%%%%%%%%%%%%%
%%% Top units %%%
%%%%%%%%%%%%%%%%%
% Top weights of fixed-N criterion (fixed-N must also be available)

% Top 10 weights of fixed-N criterion  
optimalSchemes{11}.shortName = 'top_10';
optimalSchemes{11}.longName = 'Large-N (top 10 units)';
for targetID = 1:numTargets
    optimalSchemes{11}.unrestrictedArray{targetID} = ...
        @(weightVector) boolTopCoords(weightVector, targetID, 10/numUnits);
end
optimalSchemes{11}.color =   [218, 1, 40]/255;  %  pink
optimalSchemes{11}.lineStyle = ':'; 
optimalSchemes{11}.marker = 'x';

% Top 25 weights of fixed-N criterion  
optimalSchemes{12}.shortName = 'top_25';
optimalSchemes{12}.longName = 'Large-N (top 25 units)';
for targetID = 1:numTargets
    optimalSchemes{12}.unrestrictedArray{targetID} = ...
        @(weightVector) boolTopCoords(weightVector, targetID, 10/numUnits);
end
optimalSchemes{12}.color =   [218, 1, 40]/255;  %  pink
optimalSchemes{12}.lineStyle = ':'; 
optimalSchemes{12}.marker = 'x';

% Top 10%
optimalSchemes{13}.shortName = 'top_10_pct';
optimalSchemes{13}.longName = 'Large-N (top 10\% units)';
for targetID = 1:numTargets
    optimalSchemes{13}.unrestrictedArray{targetID} = ...
        @(weightVector) boolTopCoords(weightVector, targetID, 0.1);
end
optimalSchemes{13}.color =   [90, 10, 236]/255;  % purpler
optimalSchemes{13}.colorBW = [0, 0, 0];
optimalSchemes{13}.lineStyle = '--'; 
optimalSchemes{13}.marker = 'x';
optimalSchemes{13}.markerSize = 8;

% Top 25%
optimalSchemes{14}.shortName = 'top_25_pct';
optimalSchemes{14}.longName = 'Large-N (top 25\% units)';
for targetID = 1:numTargets
    optimalSchemes{14}.unrestrictedArray{targetID} = ...
        @(weightVector) boolTopCoords(weightVector, targetID, 0.25);
end
optimalSchemes{14}.color =   [90, 10, 236]/255;  % purpler
optimalSchemes{14}.colorBW = [0, 0, 0];
optimalSchemes{14}.lineStyle = '--'; 
optimalSchemes{14}.marker = 'x';
optimalSchemes{14}.markerSize = 8;

% Top 50%
optimalSchemes{15}.shortName = 'top_50_pct';
optimalSchemes{15}.longName = 'Large-N (top 50\% units)';
for targetID = 1:numTargets
    optimalSchemes{15}.unrestrictedArray{targetID} = ...
        @(weightVector) boolTopCoords(weightVector, targetID, 0.5);
end
optimalSchemes{15}.color =   [90, 10, 236]/255;  % purpler
optimalSchemes{15}.colorBW = [0, 0, 0];
optimalSchemes{15}.lineStyle = '--'; 
optimalSchemes{15}.marker = 'x';
optimalSchemes{15}.markerSize = 8;


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Random restrictions %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Only a couple: in general random restrictions will have the same profile
% as a fixed-N regime

% Random 10 units
optimalSchemes{16}.shortName = 'random_10';
optimalSchemes{16}.longName = 'Large-N (random 10 units)';
for targetID = 1:numTargets
    randomBool = false(numUnits, 1);
    randomBool(datasample(1:numUnits, min(numUnits, 10), 'Replace',false)) = true;
    randomBool(targetID) = true;
    optimalSchemes{16}.unrestrictedArray{targetID} = ...
        @(weightVector) randomBool;
end
optimalSchemes{16}.color =   [65, 60, 174]/255;  % purplish blue
optimalSchemes{16}.lineStyle = ':'; 
optimalSchemes{16}.marker = '>';

% Random 20 units
optimalSchemes{17}.shortName = 'random_20';
optimalSchemes{17}.longName = 'Large-N (random 20 units)';
for targetID = 1:numTargets
    randomBool = false(numUnits, 1);
    randomBool(datasample(1:numUnits, min(numUnits, 20), 'Replace',false)) ...
        = true;
    randomBool(targetID) = true;
    optimalSchemes{17}.unrestrictedArray{targetID} = ...
        @(weightVector) randomBool;
end
optimalSchemes{17}.color =   [30, 54, 236]/255;  %  blue
optimalSchemes{17}.lineStyle = ':'; 
optimalSchemes{17}.marker = '<';
 

%%%%%%%%%%%%%%%%%%%%%%%
%%% Oracle cluster %%%%
%%%%%%%%%%%%%%%%%%%%%%%

optimalSchemes{18}.shortName = 'oracleClasses';
optimalSchemes{18}.longName = 'Large-N (true class labels)';
for targetID = 1:numTargets
    currentCluster = thetaLabels(targetID);
    optimalSchemes{18}.unrestrictedArray{targetID} = ...
        @(weightVector) thetaLabels == currentCluster;
end
optimalSchemes{18}.color =   [38, 185, 159]/255;  %  aquamarine
optimalSchemes{18}.lineStyle = ':'; 
optimalSchemes{18}.marker = 'v';


% Anti-oracle labels
optimalSchemes{19}.shortName = 'antiOracleClasses';
optimalSchemes{19}.longName = 'Large-N (opposite class labels)';
for targetID = 1:numTargets
    currentCluster = thetaLabels(targetID);
    optimalSchemes{19}.unrestrictedArray{targetID} = ...
        @(weightVector) thetaLabels ~= currentCluster;
end
optimalSchemes{19}.color =   [38, 100, 200]/255;  %  mildly reddish blue
optimalSchemes{19}.lineStyle = ':'; 
optimalSchemes{19}.marker = '^';

% Extract only the desired averaging schemes
optimalSchemes = optimalSchemes(averagingIncludeBool);
end



function targetBool = boolTopCoords(weightVector, targetIdx, topShare)
% boolTopCoords Returns a Boolean vector of the same length as weightVector
% The entry is true if the corresponding entry of weightVector is in the
% topShare*100 percent of values in weightVector

% Dimension of the problem
numUnits = length(weightVector);
% Extract suitable top values
numberTop = ceil(topShare*numUnits);
[~, topIds] = maxk(weightVector, numberTop);
% Create output
targetBool = false(numUnits, 1);
targetBool(topIds) = true;
targetBool(targetIdx) = true;

end