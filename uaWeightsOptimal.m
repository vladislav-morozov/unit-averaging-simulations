function optWeights = ...
    uaWeightsOptimal(estCoefs, estCovars, gradientEstimateTarget, ...
    targetID, unrestrictedBool)
% uaAveragingWeights Computes optimal weights based on estimated
% coefficients, covariances, and the supplied estimated gradient.
% The restricted units are specified in the unrestrictedBool Boolean
% vector.
%
% Args:
%     estCoefs -- kxN matrix, columns index different units
%     estCovars -- kxkxN matrix, third index indices units
%     targetIdx -- positive integer <= N, id of the target unit
%     gradientEstimateTarget -- k-vector; gradient of the target parameter
%                               with respect to thetas, evaluate at the
%                               target unit
%     unrestricted_bool: (optional) empty or Boolean vector with dimension
%                        matching estCoefs. True means the unit is
%                        unrestricted. If None, all units are unrestricted
%
% Returns:
%     optWeights -- N-vector of optimal averaging weights

% Checking for fixed-N/large-N regime
% If no value is supplied, fixed-N regime is the default
if nargin < 5 || isempty(unrestrictedBool)
    unrestrictedBool = true(size(estCoefs, 2), 1);
end

% Construct the objective function
Q = uaObjectiveMatrix(estCoefs, estCovars, gradientEstimateTarget, ...
    targetID, unrestrictedBool);
numCoords = size(Q, 1);

% Minimize
options = optimoptions('quadprog', 'Display', 'none');
xOpt = quadprog(Q, zeros(numCoords,1), [], [], ...
    ones(1, numCoords), 1, zeros(1, numCoords), [], [],...
    options);

% Fill the weights out
optWeights = nan(size(estCoefs, 2), 1); 

numRestr = length(estCoefs) - sum(unrestrictedBool);
if numRestr == 0
    % In the fixed-N regime sufficient to just spread the non-missing weights
    optWeights = xOpt;
else
    % If in large-N mode, we need to correctly allocate the weights
    % into correct non-restricted positions first, and then
    % split the residual weight across restricted units
    
    % Insert unrestricted weights
    unresWeights = xOpt(1:end-1);
    optWeights(unrestrictedBool) = unresWeights;
    
    % Split restricted weights
    resWeightInd = xOpt(end) / numRestr;
    optWeights(~unrestrictedBool) = resWeightInd;
end
end
