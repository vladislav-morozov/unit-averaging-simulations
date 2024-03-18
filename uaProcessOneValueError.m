function msePointStruct = ...
    uaProcessOneValueError(errorsArrayTarget, paramArray)
% uaProcessOneValueError Estimate MSE based on elements of the 
% errorsArrayTarget array. Uses parameter names supplied in paramArray
% 
% Returns: 
%   msePointStruct -- struct with field corresponding to parameter names.
%   Values are tables; columns are MSEs different averaging approaches

% Extract number of parameters
numParams = length(paramArray);

% Loop through parameters
for paramID = 1:numParams
    
    % Extract parameter name
    paramName = paramArray{paramID}.saveName;
    
    % Loop through sample for that parameter
    for sampleID = 1:length(errorsArrayTarget)
        
        % Initialize array if dealing with first sample
        if sampleID == 1
            sampleErrors = ...
                cell2mat(...
                struct2cell(...
                errorsArrayTarget{sampleID}.(paramName)...
                ))';
        % Otherwise add to existing array
        else
            sampleErrors = ...
                [sampleErrors; ...
                    cell2mat(...
                    struct2cell(...
                    errorsArrayTarget{sampleID}.(paramName)...
                    ))'];
        end
    end
    % Compute the MSE as the average square error
    mseEst = mean(sampleErrors.^2);
    
    % Convert to table; column names are averaging approaches
    mseEstTable = ...
        array2table(mseEst, ...
        'VariableNames', ...
        fieldnames(errorsArrayTarget{sampleID}.(paramName)));
    
    % Assign that table to a suitable output field
    msePointStruct.(paramName) = mseEstTable;
end


end