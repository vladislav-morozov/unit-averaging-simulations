function outputWeights = ...
    uaWeightsAICMMA(individualEstimators, y, x, target, k)
    % uaWeightsAICMMA Computes exponential AIC  and MMA weights to average models,
    % averaging first k models
    % target may be equal to 'aic' or 'mma'
    
    
    % Extract dimension
    [T, N,~] = size(x); 
    
    % If last argument is not supplied, use all units
    if nargin < 5
       k = N; 
    end
    
    sigmaHatSq = zeros(k,1);
    ll = zeros(k,1);
    H1 = [x(:,1,1), x(:,1,2)];
    E = zeros(T, k);
    for i=1:k
       H  = [x(:,i, 1), x(:,i,2)];
       errorsVectorI = y(1:end,i)-H*individualEstimators(:,i);
       sigmaHatSq(i) = errorsVectorI'*errorsVectorI/(T-2);
       errorsVector1 = y(1:end,1)-H1*individualEstimators(:,i);
       E(:,i) = errorsVector1;
%        ll(i) = -log(sigmaHatSq(i))/2 +1/T*(errorsVector1'*errorsVector1)/(2*sigmaHatSq(i));
       ll(i) = 1/T*(errorsVector1'*errorsVector1);
    end
    
      
    if target == "aic"
        outputWeights = exp(-ll*T + 2)/sum(exp(-ll*T + 2));
    else
        [~, minAICIdx] = max(ll);
        outputWeights = ((1:k) == minAICIdx)';
    end

    
end