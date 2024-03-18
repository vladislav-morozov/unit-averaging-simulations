function [aicW, mmaW] = uaLinearAICMMweights(individualEstimators, y, x, k)
    % aicWeights Computes exponential AIC  and MMA weights to average models,
    % averaging first k models
    
    [T, ~,~] = size(x); 
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
    
    aicW = exp(-ll*T + 2)/sum(exp(-ll*T + 2));
    options = optimoptions('quadprog','Display','none');
    mmaW = quadprog(E'*E,2*sigmaHatSq(1)*2*ones(k,1), [], [], ...
                ones(1, k), 1,zeros(1, k),[],[], options);
    
end