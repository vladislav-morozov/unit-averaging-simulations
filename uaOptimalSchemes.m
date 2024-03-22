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


% Fixed-N
optimalSchemes{1}.shortName = 'unrestr';
optimalSchemes{1}.longName = 'Fixed-N';
for targetID = 1:numTargets
    optimalSchemes{1}.unrestrictedArray{targetID} = ...
        @(weightVector) true(numUnits, 1);
end
optimalSchemes{1}.color =   [ 132, 1, 3]/255;  % dark red
optimalSchemes{1}.lineStyle = '-';
optimalSchemes{1}.marker = 'o';

% Clustering coefficients
coefClusters = kmeans(thetaHat', 2); % cluster coefs
optimalSchemes{2}.shortName = 'cluster_coef';
optimalSchemes{2}.longName = 'Large-N (clustering coefficients)';
for targetID = 1:numTargets
    currentCluster = coefClusters(targetID);
    optimalSchemes{2}.unrestrictedArray{targetID} = ...
        @(weightVector) coefClusters == currentCluster;
end
optimalSchemes{2}.color =  [90, 30, 236]/255;  % purple
optimalSchemes{2}.lineStyle = ':'; 
optimalSchemes{2}.marker = '*';

% Top weights of fixed-N criterion (fixed-N must also be available)
optimalSchemes{3}.shortName = 'top';
optimalSchemes{3}.longName = 'Large-N (top 10% unrestricted)';
for targetID = 1:numTargets
    optimalSchemes{3}.unrestrictedArray{targetID} = ...
        @(weightVector) boolTopCoords(weightVector, targetID, 0.1);
end
optimalSchemes{3}.color =   [218, 1, 136]/255;  % intense pink
optimalSchemes{3}.lineStyle = ':'; 
optimalSchemes{3}.marker = 'x';


% Stein-like
optimalSchemes{4}.shortName = 'stein';
optimalSchemes{4}.longName = 'Stein-like';
for targetID = 1:numTargets
    optimalSchemes{4}.unrestrictedArray{targetID} = ...
        @(weightVector)  ((1:numUnits)==targetID)';
end
optimalSchemes{4}.color =   [179, 112, 79]/255;  % weird brown
optimalSchemes{4}.lineStyle = ':'; 
optimalSchemes{4}.marker = 'pentagram';


% Random 10 units
optimalSchemes{5}.shortName = 'random10';
optimalSchemes{5}.longName = 'Large-N (random 10 units)';
for targetID = 1:numTargets
    randomBool = false(numUnits, 1);
    randomBool(datasample(1:numUnits, min(numUnits, 10), 'Replace',false)) = true;
    randomBool(targetID) = true;
    optimalSchemes{5}.unrestrictedArray{targetID} = ...
        @(weightVector) randomBool;
end
optimalSchemes{5}.color =   [65, 60, 174]/255;  % purplish blue
optimalSchemes{5}.lineStyle = ':'; 
optimalSchemes{5}.marker = '>';


% Random 20 units
optimalSchemes{6}.shortName = 'random20';
optimalSchemes{6}.longName = 'Large-N (random 20 units)';
for targetID = 1:numTargets
    randomBool = false(numUnits, 1);
    randomBool(datasample(1:numUnits, min(numUnits, 20), 'Replace',false)) ...
        = true;
    randomBool(targetID) = true;
    optimalSchemes{6}.unrestrictedArray{targetID} = ...
        @(weightVector) randomBool;
end
optimalSchemes{6}.color =   [30, 54, 236]/255;  %  blue
optimalSchemes{6}.lineStyle = ':'; 
optimalSchemes{6}.marker = '<';


% Random 10%
optimalSchemes{7}.shortName = 'random10pct';
optimalSchemes{7}.longName = 'Large-N (random 10%)';
for targetID = 1:numTargets
    randomBool = false(numUnits, 1);
    randomBool(...
        datasample(1:numUnits, ceil(numUnits*0.1), 'Replace',false)...
        ) = true;
    randomBool(targetID) = true;
    optimalSchemes{7}.unrestrictedArray{targetID} = ...
        @(weightVector) randomBool;
end
optimalSchemes{7}.color =   [30, 54, 236]/255;  %  blue
optimalSchemes{7}.lineStyle = ':'; 
optimalSchemes{7}.marker = '>';


% Oracle labels
optimalSchemes{8}.shortName = 'oracleClasses';
optimalSchemes{8}.longName = 'Large-N (true class labels)';
for targetID = 1:numTargets
    currentCluster = thetaLabels(targetID);
    optimalSchemes{8}.unrestrictedArray{targetID} = ...
        @(weightVector) thetaLabels == currentCluster;
end
optimalSchemes{8}.color =   [38, 185, 159]/255;  %  aquamarine
optimalSchemes{8}.lineStyle = ':'; 
optimalSchemes{8}.marker = 'v';


% Anti-oracle labels
optimalSchemes{9}.shortName = 'antiOracleClasses';
optimalSchemes{9}.longName = 'Large-N (opposite class labels)';
for targetID = 1:numTargets
    currentCluster = thetaLabels(targetID);
    optimalSchemes{9}.unrestrictedArray{targetID} = ...
        @(weightVector) thetaLabels ~= currentCluster;
end
optimalSchemes{9}.color =   [38, 100, 200]/255;  %  mildly reddish blue
optimalSchemes{9}.lineStyle = ':'; 
optimalSchemes{9}.marker = '^';


% Oracle similarity
optimalSchemes{10}.shortName = 'oracleSimilarity';
optimalSchemes{10}.longName = 'Large-N (10 most similar)';
for targetID = 1:numTargets
    dists = sum((thetaTrue-thetaTrue(:, 1)).^2);
    [~, I] = sort(dists); 
    mostSimilarIdx = I(1:10);
    similarBool = false(numUnits, 1);
    similarBool(mostSimilarIdx) = true;
    optimalSchemes{10}.unrestrictedArray{targetID} = ...
        @(weightVector) similarBool;
end
optimalSchemes{10}.color =   [0, 0, 255]/255;  %  most blue
optimalSchemes{10}.lineStyle = ':'; 
optimalSchemes{10}.marker = '<';


% Top 10 weights of fixed-N criterion  
optimalSchemes{11}.shortName = 'top10';
optimalSchemes{11}.longName = 'Large-N (top 10 unrestricted)';
for targetID = 1:numTargets
    optimalSchemes{11}.unrestrictedArray{targetID} = ...
        @(weightVector) boolTopCoords(weightVector, targetID, 10/numUnits);
end
optimalSchemes{11}.color =   [218, 1, 40]/255;  %  pink
optimalSchemes{11}.lineStyle = ':'; 
optimalSchemes{11}.marker = 'x';



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