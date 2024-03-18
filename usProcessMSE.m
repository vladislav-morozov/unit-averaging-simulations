function mseTableStruct = ... 
    usProcessMSE(msePointStructsArray, paramArray)  
% usProcessMSE Combine MSE tables in msePointStructsArray
% 
% Returns: 
%   mseTableStruct -- struct with field corresponding to parameter names.
%   Values are tables; columns are MSEs different averaging approaches

% Extract number of parameters
numParams = length(paramArray);

% Loop through parameters
for paramID = 1:numParams
    
    % Extract parameter name
    paramName = paramArray{paramID}.saveName;
    % Loop through target points
    for targetValueID = 1:length(msePointStructsArray)
        % If first value, create the table
        if targetValueID == 1
            paramMSETable = msePointStructsArray{targetValueID}.(paramName);
            % Otherwise append to existing table
        else
            paramMSETable = [...
                paramMSETable;
                msePointStructsArray{targetValueID}.(paramName)...
                ];
        end
        
        % Save to corresponding field
        mseTableStruct.(paramName) = paramMSETable;
    end
end