function [weightsTableStruct, firstUnitTable]...
    = combineAllPointsWeightsUnitsTables(paramArray, weightsRegArray, averageFirstWeightArray)

% Rows index target points, columns index grid points


% Extract number of parameters
numParams = length(paramArray);

% Loop through parameters
for paramID = 1:numParams
    
    % Extract parameter name
    paramName = paramArray{paramID}.saveName;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Processing all weights %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Extract and loop through approaches
    approachListAllWeights = fieldnames(weightsRegArray{1}.(paramName));
    % Loop through target points
    for approachID = 1:length(approachListAllWeights)
        approachName = approachListAllWeights{approachID};
        for targetValueID = 1:length(weightsRegArray)
            % If first value, create the table
            if targetValueID == 1
                weightsTable = ...
                    weightsRegArray{targetValueID}.(paramName).(approachName);
                
            else
                weightsTable = [...
                    weightsTable;
                    weightsRegArray{targetValueID}.(paramName).(approachName)...
                    ];
                
            end
            
        end
        weightsTableStruct.(paramName).(approachName) = weightsTable;
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Processing first unit weights %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Extract and loop through approaches
    approachListFirstWeight  = fieldnames(averageFirstWeightArray{1}.(paramName));
    % Loop through target points
    for approachID = 1:length(approachListFirstWeight)
        approachName = approachListFirstWeight{approachID};
        for targetValueID = 1:length(weightsRegArray)
            % If first value, create the table
            if targetValueID == 1
                firstUnitVector = ...
                    averageFirstWeightArray{targetValueID}.(paramName).(approachName);
                
            else
                firstUnitVector = [...
                    firstUnitVector;
                    averageFirstWeightArray{targetValueID}.(paramName).(approachName)...
                    ];
                
            end
            
        end
        firstUnitStruct.(paramName).(approachName) = firstUnitVector;
    end
    firstUnitTable.(paramName) = struct2table(firstUnitStruct.(paramName));
end


%%%% Convert to table only vector outputs, otherwise leave them as is


end