function thetaAvg = uaAvgEstLargeN(indEst, weights, weightsTail)
% uaAvgEstLargeN Computes a large-N regime weighted estimator thetaAvg
% using individual estimators indEst (must be a row vector) and an i0x1
% vector of weights. The tail beyond i0 is averaged with equal weights  by
% default. An alternative weight scheme can be supplied as row vector weightsTail

    i0 = length(weights);
    if nargin < 3  
        thetaAvg = indEst(1:i0-1)*weights(1:i0-1) + weights(i0)*mean(indEst(i0:end));
    else
        thetaAvg = indEst(1:i0-1)*weights(1:i0-1) + weights(i0)*(indEst(i0:end)*weightsTail);
    end
end