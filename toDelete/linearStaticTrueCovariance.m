function V = linearStaticTrueCovariance(sigmaSq, N)
    % linearStaticTrueCovariance Constructs population covariance matrices
    % for all individual estimators

    V = repmat(eye(2),1,1,N); 
    V  = reshape(sigmaSq,1,1,N).*V;
    
end