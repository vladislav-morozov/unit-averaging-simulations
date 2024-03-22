function weightsTableStruct = uaProcessWeights(weightsRegArray, paramArray)

% Rows index target points, columns index grid points


% Extract number of parameters
numParams = length(paramArray);

% Loop through parameters
for paramID = 1:numParams
    
    % Extract parameter name
    paramName = paramArray{paramID}.saveName;
    
    % Extract and loop through approaches
    approachList = fieldnames(weightsRegArray{1}.(paramName));
    % Loop through target points
    for approachID = 1:length(approachList)
        approachName = approachList{approachID};
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
end