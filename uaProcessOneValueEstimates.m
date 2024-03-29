function estPointStruct = ...
    uaProcessOneValueEstimates(estArrayTarget, paramArray)
% uaProcessOneValueEstimates Unpack estimated values for a given simulation
% target value.Uses parameter names supplied in paramArray
% 
% Returns:  

% Extract number of parameters
numParams = length(paramArray);

% Loop through parameters
for paramID = 1:numParams
    
    % Extract parameter name
    paramName = paramArray{paramID}.saveName;
    
    % Loop through sample for that parameter
    for sampleID = 1:length(estArrayTarget)
        
        % Initialize array if dealing with first sample
        if sampleID == 1
            sampleEst = ...
                cell2mat(...
                struct2cell(...
                estArrayTarget{sampleID}.(paramName)...
                ))';
        % Otherwise add to existing array
        else
            sampleEst = ...
                [sampleEst; ...
                    cell2mat(...
                    struct2cell(...
                    estArrayTarget{sampleID}.(paramName)...
                    ))'];
        end
    end
 
 
    estTable = ...
        array2table(sampleEst, ...
        'VariableNames', ...
        fieldnames(estArrayTarget{sampleID}.(paramName)));
 
    
    % Assign that table to a suitable output field
    estPointStruct.(paramName) = estTable; 
end


end