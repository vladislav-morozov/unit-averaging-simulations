function [mseTableStruct, biasTableStruct, varTableStruct] = ... 
    uaProcessMSE(msePointStructsArray,...
    biasPointStructsArray, ...
    varPointStructsArray,...
    paramArray)  
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
            paramMSETable = ...
                msePointStructsArray{targetValueID}.(paramName);
            paramBiasTable = ...
                biasPointStructsArray{targetValueID}.(paramName);
            paramVarTable = ...
                varPointStructsArray{targetValueID}.(paramName);
            % Otherwise append to existing table
        else
            paramMSETable = [...
                paramMSETable;
                msePointStructsArray{targetValueID}.(paramName)...
                ];
            paramBiasTable = [...
                paramBiasTable;
                biasPointStructsArray{targetValueID}.(paramName)...
                ];
            paramVarTable = [...
                paramVarTable;
                varPointStructsArray{targetValueID}.(paramName)...
                ];
        end
        
        % Save to corresponding field
        mseTableStruct.(paramName) = paramMSETable;
        biasTableStruct.(paramName) = paramBiasTable;
        varTableStruct.(paramName) = paramVarTable;
    end
end

end