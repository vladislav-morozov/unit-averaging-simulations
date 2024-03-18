function [thetaSample, sigmaSqSample] = ...
    uaDrawCoefficientsExp2(approach, N, varNoise, seedCoefficients)
% uaDrawCoefficients Simulates data for the AR(1) panel.
% The underlying model is
%   y_{it} = theta_{i1} y_{it-1} + theta_{i2} x_{it} + u_{it}
% Implements several options for drawing parameters (approach)
%
% Args: 1. approach -- string; the DGP for coefficients. 
%       2. N -- integer; cross-sectional size of the sample
%       3. varNoise -- integer; variance of the variances of u_{it}
%       4. seedCoefficients -- integer; seed for the rng
%
% Returns 1. thetaSamples -- 2xN matrix, columns index cross-sectional
%            units. First row holds theta_{i1}, second row holds theta_{i2}
%         2. sigmaSq -- N-vector, ith component is Var(u_{it})

% Set RNG
rng(seedCoefficients,'philox')
 
% Draw coefficients depending on the approach
if approach == "unimodal"
    % AR(1) coefficients follow a Beta distribution
    
    % Draw theta_{i1}
    lambdaMean = 0.5;
    lambdaScale = 0.7;
    lambdaSample = lambdaMean+lambdaScale*(betarnd(5, 5, 1, N)-0.5);
    
    % Draw theta_{i2}
    betaSample = 1 + sqrt(1)*randn(1, N);

elseif approach == "bimodal"
    % AR(1) coefficients follow a mixture of Betas
    % Exogeneous coefficients are N(1, 1)
    
    % Draw theta_{i1}
    lambdaSample = nan(1, N);  
    betaCenters = [0.3, 0.65];
    betaScales = [0.5, 0.3];
    for unitID = 1:N
        unitComponent = randi(2);
        lambdaSample(unitID) = betaCenters(unitComponent) +...
            betaScales(unitComponent)*(betarnd(5,5,1,1)-0.5);
    end
    
    % Draw theta_{i2}
    betaSample = 1 + sqrt(1)*randn(1, N);
    
elseif approach == "local"
    % Local parameters -- results approximately independent of T
    % Used in version 1 of the paper, all coefficients are as there
    % This setting should not be used to do comparisons across different
    % values 
    T = 60;  
    meanLambda = 0;
    meanBeta = 1;
    varianceBeta = 1;
    % Generate lambdas with locality. Support of eta scaled to compensate
    % for deflating the coefficients
    etaLambda = 2*(1-abs(meanLambda))*(2*rand([1, N])-1);  
    lambdaSample  = meanLambda+etaLambda/sqrt(T);  
    
    % Generate betas, also with locality
    betas = sqrt(varianceBeta)*randn(1, N);
    betaSample = betas/sqrt(T)+meanBeta;   
end
% Stack the two coefficient vectors together
thetaSample = [lambdaSample; betaSample];
% Draw variance of the error term
sigmaSqSample = exprnd(varNoise, 1, N);
end