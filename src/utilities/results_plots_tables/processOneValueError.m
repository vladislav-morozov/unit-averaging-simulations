function [msePointStruct, biasPointStruct, varPointStruct] = ... 
    processOneValueError(...
        errorsArrayTarget, estArrayTarget, ...
        paramArray, ...
        targetSampleParams)
% processOneValueError Estimates MSE, bias, and variance for unit averaging
% approaches.
%
% This function calculates Mean Squared Error (MSE), bias, and variance for
% each parameter in `paramArray` across multiple samples and averaging
% methods. 
%
% Args:
%     errorsArrayTarget (cell array): Cell array of structs, where each 
%         element corresponds to a sample. 
%         - Fields match the parameter names in `paramArray`. 
%         - Each field contains a struct, with keys as averaging method 
%           names and values as estimation errors.
%     estArrayTarget (cell array): Cell array of structs with the same 
%         structure as `errorsArrayTarget`. 
%         - Values are the estimated parameter values instead of errors.
%     paramArray (cell array): Array of parameter structs, where each 
%         struct includes:
%         - mu (function handle): Function to compute parameter estimates.
%         - gradient (function handle): Function to compute  gradients.
%         - shortName (string): Short identifier for the parameter.
%         - saveName (string): Key used to save results for this parameter.
%     targetSampleParams (cell array): Cell array of structs, where each 
%         element corresponds to a sample. 
%         - Fields match the parameter names in `paramArray`. 
%         - Values are true parameter values for the corresponding sample.
%
% Returns:
%     msePointStruct (struct): Struct with fields corresponding to 
%         parameter names.
%         - Each field is a table with columns representing MSEs for 
%           different averaging methods.
%     biasPointStruct (struct): Struct with the same structure as 
%         `msePointStruct`.
%         - Values in the tables are biases of the unit averaging methods.
%     varPointStruct (struct): Struct with the same structure as 
%         `msePointStruct`.
%         - Values in the tables are variances of the unit averaging
%           methods.

    % Extract number of parameters
    numParams = length(paramArray);

    % Trim 2% of outlier values when computing statistics
    meanTrimPct = 2;

    % Loop through parameters
    for paramID = 1:numParams
        % Extract parameter name
        paramName = paramArray{paramID}.saveName;

        % Initialize arrays for errors, estimates, and true values
        sampleErrors = [];
        sampleEst = [];
        sampleTarget = [];

        % Loop through samples for the current parameter
        for sampleID = 1:length(errorsArrayTarget)
            % Extract errors and estimates for the current sample
            currentErrors = ...
                cell2mat(struct2cell(errorsArrayTarget{sampleID}.(paramName)))';
            currentEst = ...
                cell2mat(struct2cell(estArrayTarget{sampleID}.(paramName)))';
            currentTarget = targetSampleParams{sampleID}.(paramName);

            % Accumulate data across samples
            sampleErrors = [sampleErrors; currentErrors];
            sampleEst = [sampleEst; currentEst];
            sampleTarget = [sampleTarget; currentTarget];
        end

        % Compute trimmed statistics
        mseEst = trimmean(sampleErrors.^2, meanTrimPct);
        biasEst = trimmean(sampleEst, meanTrimPct) - ...
                    trimmean(sampleTarget, meanTrimPct);
        varEst = var(sampleEst);

        % Convert results to tables; column names are averaging names
        colNames = fieldnames(errorsArrayTarget{1}.(paramName));
        mseEstTable = array2table(mseEst, 'VariableNames', colNames);
        biasEstTable = array2table(biasEst, 'VariableNames', colNames);
        varEstTable = array2table(varEst, 'VariableNames', colNames);

        % Assign tables to the output structs
        msePointStruct.(paramName) = mseEstTable;
        biasPointStruct.(paramName) = biasEstTable;
        varPointStruct.(paramName) = varEstTable;
    end
end
