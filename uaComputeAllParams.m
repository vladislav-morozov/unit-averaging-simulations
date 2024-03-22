function paramValues = uaComputeAllParams(paramArray, thetas, covars, y)
% uaComputeAllParams Evaluates the values of all parameters in paramArray
% using the supplied values for coefficients, covariates, and outcomes

% Number of parameters
numParams = length(paramArray);

% Loop through parameters
for paramID = 1:numParams
    % Extract parameter name
    paramName = paramArray{paramID}.saveName;
    % Compute parameter value and insert into output array
    paramValues.(paramName) = paramArray{paramID}.mu(thetas, covars, y);
end
end