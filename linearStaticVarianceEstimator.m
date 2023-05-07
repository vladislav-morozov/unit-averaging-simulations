function varianceEstimator = linearStaticVarianceEstimator(y, x, individualEstimators)
    % linearStaticTrueCovariance Constructs population covariance matrices
    % for all individual estimators


    [T, N,~] = size(x);
    varianceEstimator = zeros(2,2, N);
    sigmaHatSq = zeros(N,1);
    for i=1:N
       H = [x(:,i,1), x(:,i,2)];
       errorsVectorI = y(1:end,i)-H*individualEstimators(:,i);
       sigmaHatSq(i) = errorsVectorI'*errorsVectorI/(T-2);
       varianceEstimator(:, :, i)= (H'*H)\eye(2)*sigmaHatSq(i) ;  
    end
end