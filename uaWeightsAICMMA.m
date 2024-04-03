function outputWeights = ...
    uaWeightsAICMMA(individualEstimators, y, x, weightScheme, numUnitsAvg)
% uaWeightsAICMMA Computes exponential AIC  and MMA weights  for unit
% averaging targeting the first unit.
% AIC is computed using normal likelihood.
%
% Args: 
%       1. individualEstimators -- kxN matrix of individual estimates,
%          columns index cross-sectional units
%       2. y -- TxN matrix of outcome variables, columns index units
%       3. x -- TxNxk array of covariates of the model, third dimension
%          indexes covariates
%       4. weightScheme -- string, 'aic' or 'mma', determines which weights 
%          are returned 
%       5. k -- integer. Averaging will be done using the first k units
%
% Outputs: 
%       1. outputWeights -- k-vector of weights


% Extract dimension
[T, N, k] = size(x);

% If last argument is not supplied, use all units
if nargin < 5
    numUnitsAvg = N;
end

sigmaHatSq = zeros(numUnitsAvg,1);
ll = zeros(numUnitsAvg,1);
H1 = [x(:,1,1), x(:,1,2)]; 
for i=1:numUnitsAvg
    H  = [x(:,i, 1), x(:,i,2)];
    errorsVectorI = y(1:end,i)-H*individualEstimators(:,i);
    sigmaHatSq(i) = errorsVectorI'*errorsVectorI/(T-2);
    
    errorsVector1 = y(1:end,1)-H1*individualEstimators(:,i); 
    % Compute normal log-likelihood; note that variances are also treated
    % as part of the individual parameter vector
    ll(i) = -log(sigmaHatSq(i))/2 + ...
        1/T*(errorsVector1'*errorsVector1)/(2*sigmaHatSq(i)); 
end


if weightScheme == "aic"
    outputWeights = exp(-2*ll + 2*k)/sum(exp(-2*ll + 2*k));
else
    [~, minAICIdx] = min(ll);
    outputWeights = ((1:numUnitsAvg) == minAICIdx)';
end


end