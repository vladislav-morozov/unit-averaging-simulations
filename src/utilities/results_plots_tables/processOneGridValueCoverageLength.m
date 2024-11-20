function [coverageTable, lengthTable] = ...
    processOneGridValueCoverageLength(...
        coverageArrayTarget, lengthArrayTarget, paramArray)
% processOneGridValueCoverageLength Estimates confidence interval coverage
% and length for one grid value.
%
% This function processes Monte Carlo simulation results to compute the 
% average coverage and average length of confidence intervals for 
% various parameters and unit averaging approaches.
%
% Args:
%     coverageArrayTarget (cell array): Confidence interval coverages. 
%         - Cells correspond to Monte Carlo samples.
%         - Each cell contains a struct with the structure: 
%           paramName -> approachName -> Boolean (true if CI covers truth).
%     lengthArrayTarget (cell array): Lengths of confidence intervals. 
%         - Same structure as `coverageArrayTarget`.
%     paramArray (cell array): Definitions of parameters, where each struct
%         must include:
%         - saveName (string): Key to identify the parameter in the output.
%
% Returns:
%     coverageTable (struct): Struct of tables with estimated CI coverage.
%         - Structure: paramName -> table.
%         - Table columns correspond to unit averaging approaches.
%     lengthTable (struct): Struct of tables with estimated CI lengths.
%         - Same structure as `coverageTable`.

    % Extract the number of parameters
    numParams = length(paramArray);

    % Convert input cell arrays into structs of tables
    coverageTables = ...
        processOneSampleBootstrapEstimates(...
            coverageArrayTarget, paramArray);
    lengthTables = ...
        processOneSampleBootstrapEstimates(...
            lengthArrayTarget, paramArray);

    % Initialize output structs
    coverageTable = struct();
    lengthTable = struct();

    % Loop through each parameter
    for paramID = 1:numParams
        % Extract the parameter name
        paramName = paramArray{paramID}.saveName;

        % Estimate coverage by averaging over samples
        coverageTable.(paramName) = ...
            varfun(@mean, coverageTables.(paramName));

        % Estimate length by averaging over samples
        lengthTable.(paramName) = ...
            varfun(@mean, lengthTables.(paramName));
    end
end