function optimalWeight = uaOptimalWeights(psiMatrix)
    % averagingOptimalWeights Constructs optimal weights using asymptotic
    % covariance matrices V, individual coefficients eta, and Jacobian D
    % negativeWeights = 0, weights restricted to be nonnegative,
    %                 = 1, negative weights are allowed
    
    [~, N] = size(psiMatrix);
    
    
    if psiMatrix(1,1)~=0 % check if MG estimator has any bias
        options = optimoptions('quadprog','Display','none');
        optimalWeight = quadprog(psiMatrix,zeros(N,1), [], [], ...
                ones(1, N), 1,zeros(1, N),[],[], options);
    else
        optimalWeight = zeros(N,1); 
        optimalWeight(1) = 1;
    end
end