function [theta] = linearDynamicEstimators(y, x, N, T)
    % linearDynamicEstimators Constructs the collection of individual
    % estimators for the linear model from data 
    
    theta = zeros(2, N);
    for i=1:N
        theta(:, i) = (linearDynamicOLS(y(:,i), x(:,i), T) )';
    end
end