function fileName = ...
    makeOutputFileName(leadingStrict, simulationSetting, numReplications, valuesN, valuesT)
    % makeOutputFileName Generates a standardized file name for simulation
    % output.
    %
    % Constructs a file name string for saving simulation results based on
    % input parameters and sampler properties. The generated file name
    % follows a specific format to capture key simulation details.
    %
    % Args:
    %     thetaSampler (dataSampler): Instance of the dataSampler class 
    %         representing the distribution for theta coefficients.
    %     uSampler (dataSampler): Instance of the dataSampler class 
    %         representing the distribution for u (error terms).
    %     N (int): Cross-sectional sample size (number of individuals).
    %     T (int): Number of time periods (cross-sections).
    %     numSamples (int): Number of Monte Carlo samples in the simulation
    %     simContext (str): String identifier for the simulation context 
    %         (e.g., 'baseline' or other experiment descriptors).
    %
    % Returns:
    %     fileName (str): Constructed file name including all relevant 
    %         parameters, saved in the 'outputs' folder, with format:
    %         'outputs/[simContext]_samples_[numSamples]_N_[N]_T_[T]_...
    %          F_[thetaDistr]_[thetaParam]_[thetaParamValue]_...
    %          G_[uDistr]_[uParam]_[uParamValue].mat'

    
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
        leadingStrict, ...          % Insert a leading string
        simulationSetting, ...      % Simulation context
        numReplications,    ...     % Number of samples
        titleN, titleT);            % Sample sizes and cross-sections 
end