function V = uaTrueAsymptoticVariance(lambda, beta, sigmaSq, N, varianceX)
    % uaTrueAsymptoticVariance Constructs asymptotic covariance matrices
    % for all individuals. This is the usual \sigma^2\E(X'X)^{-1}

    % V is the covariance matrix, we immediately invert x_1x_1' due to
    % simplicity
    V = repmat(eye(2)/varianceX,1,1,N);
    V(1,1,:) = (1-reshape(lambda.^2,1,1,N))./(reshape(beta.^2,1,1,N)*varianceX + reshape(sigmaSq,1,1,N))   ;
    V  = reshape(sigmaSq,1,1,N).*V;
    
end