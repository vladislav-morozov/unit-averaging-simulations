function [betaMatrix, estVarianceArray, estVarianceUVector] = ...
    uaOLS(y, x, estimationWindow)
% uaOLS Regress y on x respecting the data structure of the simulation.
% Also estimates the variance-covariance matrix of the estimates
%
% Inputs: 1. y -- NxT matrix. Columns index individual cross-sectional
%                 units. Rows index observations.
%         2. x -- NxTxk matrix. Third dimension indexes covariates
%         3. estimationWindow -- (optional) positive integer. If specified,
%                                only the last estimationWindow
%                                observations will be used for estimation
%
% Outputs: 1. betaMatrix -- kxN matrix of OLS estimates. jth column
%                           contains the OLS regression coefficients of the
%                           jth column of y on the jth columns of x
%          2. estVarianceArray -- kxkxN array, where the third axis indexes
%                                 cross-sectional units. The (j, j)th slice
%                                 across the third dimension are the
%                                 estimated covariance matrix of
%                                 betaMatrix(:, j)
%          3. estVarianceUVector -- N-vector. jth component estimates the
%                                   variance of u_{jt}

% Use all data if estimationWindow is not supplied. 
if nargin<3
   estimationWindow = size(y, 1);
end

% Extract dimensions
[~, N, k] = size(x);

% Allocate array for coefficient estimates
betaMatrix = nan(k, N);
% Allocate variance array
estVarianceArray = nan(k, k, N);
% Allocate auxiliary array for variances of u
estVarianceUVector = nan(N,1);


% Compute estimates for each unit individually
for i=1:N
    % Extract individual data
    yInd = y(end+1-estimationWindow:end, i);
    xInd = squeeze(x(end+1-estimationWindow:end, i, :));
    
    % Compute coefficient estimates
    betaMatrix(:, i) = xInd\yInd;
    
    % Estimate variance
    errorsVectorI = yInd-xInd*betaMatrix(:,i);
    estVarianceUVector(i) = ...
        errorsVectorI'*errorsVectorI/(estimationWindow-2);
    estVarianceArray(:, :, i)= (xInd'*xInd)\eye(k)*estVarianceUVector(i) ;  
end

end