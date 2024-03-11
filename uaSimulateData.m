function [y, covars, u] = ...
    uaSimulateData(thetasMatrix, varianceUVector, varX,  ...
                    T, seedData)
% uaSimulateData Simulates data from an AR(1) model with one
% exogenous covariate. The model is:
%       y_{it} = theta_{i1} y_{it-1} + theta_{i2}x_{it} + u_{it}
%
% Draws and returns a dataset using specified coefficients theta.
% The number of cross-sectional units is determined by the second
% dimension of the thetasMatrix
%
% Inputs:
% 1. thetasMatrix -- 2xN matrix. The ith column is of the form
%                    (theta_{i1}, theta_{i2}).
% 2. varianceUVector -- N-vector of individual error variances. ith
%                       component is the variance of u_{it}
% 3. varX -- variance of the exogenous covariate x
% 4. T -- number of time periods to draw
% 5. seedData -- seed for rng
%
% Outputs
% 1. y -- TxN matrix of outcomes. Columns correspond to different units
% 2. covars -- observed covariates. x is TxNx2, first slice houses lag of y,
%              second houses exogenous variable x
% 3. u -- TxN matrix of error terms.


% Set data seed according to the argument
rng(seedData,'philox')

% Detect implied cross-section size
[~, N] = size(thetasMatrix);

% Allocate vector
y = nan(T, N);
covars = nan(T, N, 2);

% The first slice of covars across the third dim is the lag of y
% We first fill it out with the initial value that ensures covariance
% stationarity of the process
covars(1, :, 1) = mvnrnd(zeros(N, 1), ...
    diag((varX + varianceUVector)./(1-thetasMatrix(1,:).^2)),...
    1); 

% The second slice of covars is the exogeneous covariate x, iid normal over
% time with variance varX
covars(:, :, 2) = mvnrnd(zeros(N, 1), varX*eye(N), T);

% u is multivariarte normal with unit-specified variances specified by 
% varianceUVector
u = mvnrnd(zeros(N, 1), diag(varianceUVector), T);

% Simulate the AR(1) process for y
for t=1:T
    % Fill out the outcome variable
    y(t,:) = thetasMatrix(1,:).*covars(t,:,1) + ...
             thetasMatrix(2,:).*covars(t,:, 2) + ...
             u(t,:);
         
    % Fill out the lag of the variable
    if t<T
        covars(t+1, :, 1) = y(t, :);
    end
end
end