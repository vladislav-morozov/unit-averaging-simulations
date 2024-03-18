function [y, x, u] = linearDynamicSimulateData(eta, sigmaSq, T, design, seedData)
    % linearDynamicSimulateData Simulates data from an AR(1) with one
    % exogenous covariate
    % The model is y_{it} = \lambda_{i} y_{it-1} + \beta_{i}x_{it} + u_{it}
    % Draws data according to design chosen and returns data and parameter
    % vectors. Designs values accepted -- '1', '2'
    
    [~, N] = size(eta);
    y = zeros(T, N);
%     y(1,:) = y0; % put the initial conditions in y
    x = zeros(T, N);

    if design==1 % Draw data
        rng(seedData,'philox')
        x(1, :, 1) = mvnrnd(zeros(1, 1), 1, N); 
        x(:, :, 2) = mvnrnd(ones(1, N), eye(N), T); 
        u = mvnrnd(zeros(N, 1), (diag(sigmaSq)), T);
    else 
        rng(seedData,'philox')
        x(1, :, 1) = mvnrnd(eta(2,:).^2./(1-eta(1,:)),eye(N),1); 
        x(:, :, 2) = mvnrnd(eta(2,:), eye(N), T); 
        u = mvnrnd(zeros(N, 1), (diag(sigmaSq)), T);     
    end

    % Simulate data
    for t=1:T
        y(t,:) = eta(1,:).*x(t,:,1) + eta(2,:).*x(t,:, 2) + u(t,:) ;
        if t<T
        x(t+1, :, 1) = y(t, :);
        end
    end
end