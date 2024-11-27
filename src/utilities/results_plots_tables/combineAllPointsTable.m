function [combinedTableStruct] = ...
    combineAllPointsTable(sampleStructsArray, paramArray)
% combineAllPointsTable Combines tables across all grid points of theta
% into a single table.
%
% This function merges tables corresponding to different grid points for
% each parameter into a single comprehensive table. Can be used for MSE,
% bias, variance, coverage, and interval lengths.
%
% Args:
%     sampleStructsArray (cell array): Each cell corresponds to a grid 
%         point value for theta. Each cell contains a struct where:
%         - Keys are parameter names.
%         - Values are tables of results for those parameters. Columns 
%           represent different unit averaging approaches.
%     paramArray (cell array): Array of structs defining parameters. Each
%          struct must include:  
%         - `saveName` (string): The key used to identify parameter results
%
% Returns:
%     combinedTableStruct (struct): Struct where:
%         - Keys are parameter names.
%         - Values are combined tables across all grid points for the 
%           corresponding parameter.
%         - Each table's rows represent grid point values for theta, and
%            columns correspond to different unit averaging approaches.

% Extract the total number of parameters to process.
numParams = length(paramArray);

% Loop through each parameter to aggregate data.
for paramID = 1:numParams
    % Retrieve the name of the parameter from the parameter array.
    paramName = paramArray{paramID}.saveName;

    % Initialize an empty table for the current parameter.
    paramTable = [];

    % Loop through all grid points to gather data for the current parameter.
    for targetValueID = 1:length(sampleStructsArray)
        % Retrieve the table for the current grid point and parameter.
        currentTable = sampleStructsArray{targetValueID}.(paramName);

        % Append the table for the current grid point to the cumulative
        % table.
        paramTable = [paramTable; currentTable]; %#ok<AGROW>
    end

    % Store the combined table in the output struct using the parameter
    % name as the key.
    combinedTableStruct.(paramName) = paramTable;
end
end
