function [betaMatrix, estVarianceArray, estVarianceUVector] = ...
    OLS(y, x, estimationWindow)
% OLS Performs Ordinary Least Squares (OLS) regression of y on x and
% estimates the variance-covariance matrix of the estimates.
%
% Args:
%     y (matrix): An NxT matrix where columns index individual 
%         cross-sectional units and rows index observations.
%     x (3D matrix): An NxTxk matrix where the third dimension indexes 
%         covariates.
%     estimationWindow (optional, int): A positive integer specifying the 
%         number of last observations to be used for estimation. If not 
%         specified, all observations are used.
%
% Returns:
%     betaMatrix (matrix): A kxN matrix of OLS estimates. The jth column
%         contains the OLS regression coefficients of the jth column of y 
%         on the jth block of x.
%     estVarianceArray (cell array): An N-cell array where each cell 
%         contains the estimated covariance matrix of betaMatrix(:, j).
%     estVarianceUVector (vector): An N-vector where each component 
%         estimates thevariance of u_{jt}.
%
% Example:
%     y = randn(100, 5);  % 100 observations for 5 units
%     x = randn(100, 5, 3);  % 3 covariates
%     estimationWindow = 50;
%     [betaMatrix, estVarianceArray, estVarianceUVector] = ...
%       OLS(y, x, estimationWindow);
%
% Note:
%     If estimationWindow is not supplied, all available observations are
%     used for estimation.

% Use all data if estimationWindow is not supplied
if nargin < 3
    estimationWindow = size(y, 1);
end

% Extract dimensions
[~, N, k] = size(x);

% Allocate array for coefficient estimates
betaMatrix = nan(k, N);
% Allocate variance array
estVarianceArray = cell(N, 1);
% Allocate auxiliary array for variances of u
estVarianceUVector = nan(N, 1);

% Compute estimates for each unit individually
for i = 1:N
    % Extract individual data
    yInd = y(end+1-estimationWindow:end, i);
    xInd = squeeze(x(end+1-estimationWindow:end, i, :));
    
    % Compute coefficient estimates
    betaMatrix(:, i) = xInd \ yInd;
    
    % Estimate variance
    errorsVectorI = yInd - xInd * betaMatrix(:, i);
    estVarianceUVector(i) = errorsVectorI' * errorsVectorI / (estimationWindow - 2);
    estVarianceArray{i} = (xInd' * xInd) \ eye(k) * estVarianceUVector(i);
end
end