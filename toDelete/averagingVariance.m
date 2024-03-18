function AMSE = averagingVariance(weights, Psi)
    % averagingVariance Computes variance of the averaging estimator. Needs
    % weights, asymptotic covariance of the estimators, and the Jacobian D
    % of the true parameter evaluated at unit 1
    
    N = length(weights)-1;
    weights = reshape(weights, N+1,1);
    AMSE = weights'*Psi*weights;
end