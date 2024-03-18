function [beta] = linearStaticEstimators(y, x)
    % linearDynamicEstimators Constructs the collection of individual
    % estimators for the linear model from data 
    
    [~, N,~] = size(x);
    beta = zeros(2, N);
    for i=1:N
        beta(:, i) = (linearStaticOLS(y(:,i), [x(:,i,1), x(:,i,2)]) );
    end
end