function V = linearDynamicTrueCovariance(lambda, beta, sigmaSq, N, design)
    % linearDynamicTrueCovariance Constructs population covariance matrices
    % for all individual estimators

    V = repmat(eye(2),1,1,N);
    if design==1
       V(1,1,:) = (1-reshape(lambda.^2,1,1,N))./(reshape(beta.^2,1,1,N) + reshape(sigmaSq,1,1,N))   ;
    elseif design == 2
       V(1,1,:) = (1-reshape(lambda.^2,1,1,N))./(reshape(beta.^2,1,1,N).*(1+reshape(beta.^2,1,1,N)) + reshape(sigmaSq,1,1,N))   ;
    end
    
    V  = reshape(sigmaSq,1,1,N).*V;
    
end