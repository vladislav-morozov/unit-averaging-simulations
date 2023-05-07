function [y, x, u] = uaSimulateData(theta, sigmaSq,varX,  T, seedData)
    % uaSimulateData Simulates data from an AR(1) model with one
    % exogenous covariate
    % The model is y_{it} = \lambda_{i} y_{it-1} + \beta_{i}x_{it} + u_{it}
    % Draws and returns a dataset.
    %
    % Inputs: 
    % 1. theta -- individual parameter vectors, each column (of N) a unit
    % 2. sigmaSq -- N-vector of individual error variances
    % 3. varX -- variance of X
    % 4. T -- number of time periods to draw
    % 5. seedData -- seed for rng
    %
    % Outputs
    % 1. y, x -- data
    % 2. u -- error terms
    
    rng(seedData,'philox')
    [~, N] = size(theta);
    % Create vectors
    y = zeros(T, N);
    x = zeros(T, N);

    % Coordinate 1 is the lag of y
%     x(1, :, 1) = mvnrnd(theta(2,:)./(1-theta(1,:)),diag((varX + sigmaSq)./(1-theta(1,:).^2)),1); % Draw initial condition
    x(1, :, 1) = mvnrnd(zeros(N, 1),diag((varX + sigmaSq)./(1-theta(1,:).^2)),1); % Draw initial condition    
    x(:, :, 2) = mvnrnd(zeros(N, 1), varX*eye(N), T);
    u = mvnrnd(zeros(N, 1), diag(sigmaSq), T);
     
    % Simulate data
    for t=1:T
        y(t,:) = theta(1,:).*x(t,:,1) + theta(2,:).*x(t,:, 2) + u(t,:) ;
        if t<T
        x(t+1, :, 1) = y(t, :);
        end
    end
end