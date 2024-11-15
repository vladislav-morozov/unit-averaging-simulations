function paramValues = computeAllParams(paramArray, thetas, covars, y)
% computeAllParams Computes the values of all parameters in paramArray
% using the supplied values for coefficients (thetas), covariates (covars), 
% and outcomes (y).
%
% Args:
%     paramArray (cell array): A cell array of structs. Each element 
%         defines a parameter and must contain a field called 'mu', which 
%         is a functionhandle that computes the parameter value. The 
%         function handle must accept thetas, covars, and y as inputs.
%     thetas (vector): A vector of parameter values. Each element in thetas 
%         must also have a field called 'paramName' containing the name of 
%         the parameter.
%     covars (matrix): A matrix of covariate values.
%     y (vector): A vector of outcome values.
%
% Returns:
%     paramValues (struct): A struct where each field is named according to 
%         'paramName' in paramArray and contains the corresponding computed 
%          parameter value.
%
% Example:
%     numParams = 2;
%     paramArray = {struct('mu', @(theta, covar, y) mean(y), ...
%                    'paramName', 'yMean')};
%     thetas = [0.5, 1.2];
%     covars = rand(100, 2);
%     y = rand(100, 1);
%     paramValues = computeAllParams(paramArray, thetas, covars, y);
%
% Note:
%     The requirements for thetas, covars, and y are determined by the 
%     requirements of the functions in the 'mu' field of paramArray.

% Number of parameters
numParams = length(paramArray);

% Initialize output struct
paramValues = struct();

% Loop through each parameter
for paramID = 1:numParams
    % Extract the name of the current parameter
    paramName = paramArray{paramID}.paramName;
    
    % Compute the parameter value using the specified function 'mu'
    paramValues.(paramName) = paramArray{paramID}.mu(thetas, covars, y);
end
end