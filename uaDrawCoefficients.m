function [theta, thetaScaled, sigmaSq] = ...
            uaDrawCoefficients(N, T, meanCoef, varianceBeta, varNoise, seedCoefficients)
    % linearDynamicSimulateData Simulates data for the first MC experiment
    % The model is y_{it} = \lambda_{i} y_{it-1} + \beta_{i}x_{it} + u_{it}
    % Draws parameters with lambda uniform symmetrically around its mean,
    % \beta normal; variance of u_{it} exponential with parameter lambda
    
    
    
    % Draw coefficients independently
    rng(seedCoefficients,'philox')
    
    meanLambda = meanCoef(1);
    meanBeta = meanCoef(2);
    
    % Generate lambdas
    etaLambda = (1-abs(meanLambda))*(2*rand([1, N])-1);  % automatically ensure that support within [-1, 1]
    lambda = meanLambda + etaLambda; 
    lambdaScaled  = meanLambda+etaLambda/sqrt(T); % scale according to locality
    
    % Generate betas
    betas = sqrt(varianceBeta)*randn(1, N);
    betaScaled = betas/sqrt(T)+meanBeta; % scale according to locality
    betas = betas  + meanBeta;
    
    % Draw variance of the error term
    sigmaSq = exprnd(varNoise, 1, N);
    
    
    theta  = [lambda; betas];
    thetaScaled = [lambdaScaled; betaScaled];
    
    

end