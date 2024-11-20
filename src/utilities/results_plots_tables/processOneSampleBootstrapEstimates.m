function bsEstTableStruct = ...
    processOneSampleBootstrapEstimates(estArrayTarget, paramArray)
% processOneSampleBootstrapEstimates Converts bootstrap parameter
% estimates into structured tables.
%
% This function processes bootstrap estimates for multiple parameters and
% averaging approaches, organizing the results into a structured format for
% further analysis.
%
% Args:
%     estArrayTarget (cell array): Bootstrap estimates for all parameters 
%         across all samples. Each element is a struct of the form: 
%         paramName -> averagingApproach -> debiased estimate.
%     paramArray (cell array): Array of parameter definitions, where each 
%         struct includes:
%         - mu (function handle): Function to compute parameter estimates.
%         - gradient (function handle): Function to compute gradients.
%         - shortName (string): Short identifier for the parameter.
%         - saveName (string): Key used to save results for this parameter.
%
% Returns:
%     bsEstTableStruct (struct): Struct containing tables of debiased unit
%         averaging estimates.
%         - Structure: paramName -> table of debiased estimates.
%         - Table columns correspond to different unit averaging approaches

    % Extract the number of parameters
    numParams = length(paramArray);

    % Initialize output struct
    bsEstTableStruct = struct();

    % Loop through each parameter
    for paramID = 1:numParams
        % Extract parameter name from the paramArray
        paramName = paramArray{paramID}.saveName;

        % Initialize an array to store estimates across samples
        sampleEst = [];

        % Loop through each bootstrap sample
        for sampleID = 1:length(estArrayTarget)
            % Extract estimates for the current parameter and sample
            sampleData = struct2cell(estArrayTarget{sampleID}.(paramName));
            sampleData = cell2mat(sampleData)'; % Convert to row vector

            % Append the estimates to the growing array
            if sampleID == 1
                % First sample initializes the array
                sampleEst = sampleData;
            else
                % Concatenate subsequent samples
                sampleEst = [sampleEst; sampleData]; 
            end
        end

        % Convert the aggregated array of estimates into a table
        estTable = array2table(sampleEst, ...
            'VariableNames', fieldnames(estArrayTarget{1}.(paramName)));

        % Assign the table to the corresponding field in the output struct
        bsEstTableStruct.(paramName) = estTable;
    end
end
