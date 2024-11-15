function [y, covars, u] = ...
    drawData(thetasMatrix, varianceUVector, varX, T, seedData)
    % drawData Simulates data from an AR(1) model with one exogenous
    % covariate.
    %
    % This function generates data based on an AR(1) process with one
    % exogenous covariate, governed by the model:
    %       y_{it} = theta_{i1} * y_{it-1} + theta_{i2} * x_{it} + u_{it}
    %
    % Args:
    %     thetasMatrix (matrix): 2 x N matrix, where each column represents
    %         the coefficients for the AR(1) model of a particular unit
    %         (theta_{i1}, theta_{i2}).
    %     varianceUVector (vector): N-element vector, where each component 
    %         represents the variance of the error term for each unit.
    %     varX (scalar): Variance of the exogenous covariate x, assumed to 
    %         be the same for all units.
    %     T (int): Number of time periods to generate data for.
    %     seedData (int): Seed for the random number generator to ensure 
    %         reproducibility.
    %
    % Returns:
    %     y (matrix): T x N matrix of outcomes. Each column represents the 
    %         outcome series for one unit.
    %     covars (3D array): T x N x 2 array of observed covariates. 
    %         The first slice (covars(:,:,1)) holds the lagged y values,
    %         and the second slice (covars(:,:,2)) holds the exogenous 
    %         covariate x.
    %     u (matrix): T x N matrix of error terms, where each column 
    %         corresponds to a unit.
    %
    % Example:
    %     [y, covars, u] = ...
    %       drawData(thetasMatrix, varianceUVector, varX, 100, 42);

    % Set random seed to ensure reproducibility
    rng(seedData, 'philox')

    % Determine the number of cross-sectional units (N) from thetasMatrix
    [~, N] = size(thetasMatrix);

    % Initialize output matrices with NaNs
    y = nan(T, N);
    covars = nan(T, N, 2);

    % Generate the initial lagged values of y to ensure covariance
    % stationarity First slice (covars(:,:,1)) is the lagged value of y;
    % initialize with stationary values
    covars(1, :, 1) = mvnrnd(zeros(N, 1), ...
        diag((varX + varianceUVector) ./ (1 - thetasMatrix(1, :).^2)), 1);

    % Generate the exogenous covariate x (second slice of covars,
    % covars(:,:,2)) Assumed to be iid normal over time with variance varX
    covars(:, :, 2) = mvnrnd(zeros(N, 1), varX * eye(N), T);

    % Generate the error terms u as a multivariate normal with
    % unit-specific variances
    u = mvnrnd(zeros(N, 1), diag(varianceUVector), T);

    % Simulate the AR(1) process for y using the lagged values and the
    % exogenous covariate x
    for t = 1:T
        % Calculate the outcome y for each unit based on current covariate
        % values and errors
        y(t, :) = thetasMatrix(1, :) .* covars(t, :, 1) + ...    % lag y 
                  thetasMatrix(2, :) .* covars(t, :, 2) + ...    % x
                  u(t, :);                                       % Error 

        % Update lagged y values for next time period (if within bounds)
        if t < T
            covars(t + 1, :, 1) = y(t, :);
        end
    end
end
