function [beta, betaScaled, sigmaSq] = ...
            linearStaticDrawCoefficients(N, T, meanBeta, varianceBeta,...
            correlationBeta, seedCoefficients)

   
    

        % Draw coefficients independently
        rng(seedCoefficients,'philox')
        C = [varianceBeta, correlationBeta*varianceBeta; correlationBeta*varianceBeta, varianceBeta];
        beta = mvnrnd(meanBeta, C, N)';
        % Scale them locally
        betaScaled = (beta-repmat(meanBeta,1,N))/sqrt(T)+repmat(meanBeta,1,N);
        % Draw variance of the error term
        rng(seedCoefficients,'philox')
        sigmaSq = exprnd(1, 1, N);
        

   
end