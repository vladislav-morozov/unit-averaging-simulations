function optimalSchemes = uaOptimalSchemes(thetaHat)
% uaOptimalSchemes Implements the fixed-N and various large-N optimal
% schemes. 
%
% Inputs: thetaHat -- kxN matrix of coefficient estimates. Columns index
%   units
%
% Returns: optimalSchemes -- cell array of structs. Each struct describes a
%   an averaging approaches. An approach is characterized by 3 fields:
%   .shortName, .longName, and .unrestrictedArray. The latter is an array 
%   of functions. The jth function describes the way to determine the
%   unrestricted units for averaging the jth unit with a given scheme.


numUnits = size(thetaHat, 2);

% Do any sample-specific large-N computations
optimalSchemes = {};


% Fixed-N
optimalSchemes{1}.shortName = 'unrestr';
optimalSchemes{1}.longName = 'Fixed-N';
for targetID = 1:numUnits
    optimalSchemes{1}.unrestrictedArray{targetID} = ...
        @(weightVector) true(numUnits, 1);
end


% Random units
optimalSchemes{2}.shortName = 'random';
optimalSchemes{2}.longName = 'Large-N (random 10%)';
for targetID = 1:numUnits
    randomBool = false(numUnits, 1);
    randomBool(datasample(1:numUnits, ceil(0.1*numUnits))) = true;
    randomBool(targetID) = true;
    optimalSchemes{2}.unrestrictedArray{targetID} = ...
        @(weightVector) randomBool;
end


% Clustering coefficients
coefClusters = kmeans(thetaHat', 4); % cluster coefs

optimalSchemes{3}.shortName = 'cluster_coef';
optimalSchemes{3}.longName = 'Large-N (clustering coefficients)';
for targetID = 1:numUnits
    currentCluster = coefClusters(targetID);
    optimalSchemes{3}.unrestrictedArray{targetID} = ...
        @(weightVector) coefClusters == currentCluster;
end


% Top weights of fixed-N criterion (fixed-N must also be available)
optimalSchemes{4}.shortName = 'top';
optimalSchemes{4}.longName = 'Large-N (top 10% unrestricted)';
for targetID = 1:numUnits
    optimalSchemes{4}.unrestrictedArray{targetID} = ...
        @(weightVector) boolTopCoords(weightVector, targetID, 0.1);
end


% Stein-like
optimalSchemes{5}.shortName = 'stein';
optimalSchemes{5}.longName = 'Stein-like';
for targetID = 1:numUnits
    optimalSchemes{5}.unrestrictedArray{targetID} = ...
        @(weightVector)  (1:numUnits)==targetID;
end

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