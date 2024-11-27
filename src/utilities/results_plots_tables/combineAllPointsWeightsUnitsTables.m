function [combinedStruct]...
    = combineAllPointsWeightsUnitsTables(sampleArray, paramArray)
% combineAllPointsWeightsUnitsTables Combines non-table results into
% a coherent object
%
% This function merges arrays corresponding to different grid points for
% each parameter into a single comprehensive array.
% Can be used for weight vectors, average own weights, mass differences
% between restricted and unrestricted units, and probabilities of being
% unrestricted
%
% Args:
%     sampleStructsArray (cell array): Each cell corresponds to a grid 
%         point value for theta. Each cell contains a struct where:
%         - Keys are parameter names.
%         - Values are structs indexed by unit averaging approaches. Values
%           of those fields --- vectors of results 
%     paramArray (cell array): Array of structs defining parameters. Each
%          struct must include:  
%         - `saveName` (string): The key used to identify parameter results
%
% Returns:
%     combinedTableStruct (struct): Struct where:
%         - Keys are parameter names.
%         - Values are are combined results. 
%         - Results returned as a table if the individual result for
%           a given grid point and averaging approach is scalar
%         - Results returned as a struct indexed by averaging approaches if
%           the individual result for a given grid point and averaging
%           approachs is a vector

% Extract number of parameters
numParams = length(paramArray);

% Loop through parameters
for paramID = 1:numParams
    
    % Extract parameter name
    paramName = paramArray{paramID}.saveName;
       
    % Extract and loop through approaches
    approachListAllWeights = fieldnames(sampleArray{1}.(paramName));

    % Loop through approaches
    for approachID = 1:length(approachListAllWeights)

        % Extract approach name
        approachName = approachListAllWeights{approachID};
        
        % Create temporary results array
        tempResultsArray = [];

        % Loop through grid points for theta
        for targetValueID = 1:length(sampleArray)

            tempResultsArray = [...
                    tempResultsArray;
                    sampleArray{targetValueID}.(paramName).(approachName)...
                    ];
            
        end
        
        % Insert results into the output struct
        combinedStruct.(paramName).(approachName) = tempResultsArray;
    end

    % If results are scalars, convert output to a table. Otherwise do not 
    if any(size(tempResultsArray)==1)
        combinedStruct.(paramName) = ...
            struct2table(combinedStruct.(paramName));
    end

end
end