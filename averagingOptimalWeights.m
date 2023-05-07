function optimalWeight = averagingOptimalWeights(psiMatrix, negativeWeights)
    % averagingOptimalWeights Constructs optimal weights using asymptotic
    % covariance matrices V, individual coefficients eta, and Jacobian D
    % negativeWeights = 0, weights restricted to be nonnegative,
    %                 = 1, negative weights are allowed
    
    [~, N] = size(psiMatrix);
   
    assert(ismember(negativeWeights,[1,0]),'Invalid weight option')
    
    
    if psiMatrix(1,1)~=0 % check if MG estimator has any bias
        options = optimoptions('quadprog','Display','none');
        if negativeWeights == 0 % Restrict weights to be nonnegative
            optimalWeight = quadprog(psiMatrix,zeros(N,1), [], [], ...
                ones(1, N), 1,zeros(1, N),[],[], options);
        else 
            optimalWeight = quadprog(psiMatrix,zeros(N,1), [], [], ...
                ones(1, N), 1,[],[],[], options);    
        end
    else
        optimalWeight = zeros(N,1); 
        optimalWeight(1) = 1;
    end
end