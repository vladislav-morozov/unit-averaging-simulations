function [coverageTable, lengthTable] = ...
    uaProcessCoverageLength(coverageArrayTarget, lengthArrayTarget, paramArray)

% Combine iterations into a table
coverageTables = ...
    uaProcessOneValueEstimates(coverageArrayTarget, paramArray);
lengthTables = ...
    uaProcessOneValueEstimates(lengthArrayTarget, paramArray);

numParams = length(paramArray);

% Loop through parameters
for paramID = 1:numParams
    % Extract parameter name
    paramName = paramArray{paramID}.saveName;
    
    % Compute coverage
    coverageTable.(paramName) = varfun(@mean, coverageTables.(paramName));
    % Compute average length
    lengthTable.(paramName) = varfun(@mean, lengthTables.(paramName));
    
end

end