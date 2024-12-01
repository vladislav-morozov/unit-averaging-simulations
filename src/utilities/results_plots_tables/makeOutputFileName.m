function fileName = ...
    makeOutputFileName(leadingString, simulationSetting, numReplications, ...
                       valuesN, valuesT)
    % makeOutputFileName Generates a standardized file name for simulation
    % output.
    %
    % Constructs a file name string for saving simulation results based on
    % input parameters and sampler properties. The generated file name
    % follows a specific format to capture key simulation details.
    %
    % Args:
    %     leadingString (string): String to insert at the beginning of the
    %         file name.
    %     simulationSetting (string): Name of the data-generating process.
    %     numReplications (int): Number of Monte Carlo samples drawn.
    %     valuesN (vector): Vector of cross-sectional sizes used.
    %     valuesT (vector): Vector of time series sizes used.
    %
    % Returns:
    %     fileName (str): Constructed file name including all relevant 
    %         parameters, saved in the 'results/simulation' folder. Format:
    %         'results/simulation/[leadingString][simulationSetting]_...
    %         [numReplications]_N-[valuesN]_T-[valuesT].mat'
    
    % Numbers of dimensions of samples drawn
    numT = length(valuesT);
    numN = length(valuesN);

    % Generate parts of the title that depend on (N, T)
    titleN = "";
    for nID=1:numN
        titleN = titleN + "-" + valuesN(nID) ;
    end

    titleT = "";
    for tID=1:numT
        titleT = titleT + "-" + valuesT(tID);
    end

    % Construct the file name using strings with specified structure
    fileName = sprintf('results/simulations/%s%s_%d_N%s_T%s.mat', ...
        leadingString, ...          % Insert a leading string
        simulationSetting, ...      % Simulation context
        numReplications,    ...     % Number of samples
        titleN, titleT);            % Sample sizes and cross-sections 
end