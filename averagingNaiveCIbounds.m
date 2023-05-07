function [bias, std] = averagingNaiveCIbounds(thetaHat, weights, D, V, T)
    % averagingNaiveCIbounds Computes the native bound treating weights as
    % fixed
    
    variances = 0;
    variancesBias =0;
    bias = 0;
    V = V*T;
    
    for i=1:length(weights)
        variances = variances+  weights(i)^2*D'*V(:,:,i)*D; % self
        if i>1
          variancesBias = variancesBias+weights(i)^2*D'*(V(:,:,i)+V(:,:,1))*D;
        end 
        bias = bias + weights(i)*D'*sqrt(T)*(thetaHat(:,i)-thetaHat(:,1));
    end
    

    bias= bias/sqrt(T);
    std = 1.96*sqrt(variances/T)+7*sqrt(variancesBias/T);
    
end