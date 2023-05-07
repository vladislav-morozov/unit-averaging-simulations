function [theta, thetaScaled, sigmaSq] = ...
            linearDynamicDrawCoefficients(N, T, design, meanBeta, varianceBeta, seedCoefficients)
    % linearDynamicSimulateData Simulates data for the first MC experiment
    % The model is y_{it} = \lambda_{i} y_{it-1} + \beta_{i}x_{it} + u_{it}
    % Draws data according to design chosen and returns data and parameter
    % vectors. Designs values accepted -- '1', '2'
    
    
    assert(ismember(design,[1,2]),'Invalid design choice')
    
    
    
    % Draw data and parameters

    if design<=2
        % Draw coefficients independently
        rng(seedCoefficients,'philox')
%         lambda = rand([1,N])/2+.25;
        lambda = 2*rand([1, N])-1;
        rng(seedCoefficients,'philox')
        beta = sqrt(varianceBeta)*randn(1, N)+meanBeta;
        % Scale them locally
        lambdaScaled = (lambda-.5)/sqrt(T)+0.5;
        betaScaled = (beta-meanBeta)/sqrt(T)+meanBeta;
        % Draw variance of the error term
        rng(seedCoefficients,'philox')
        sigmaSq = exprnd(1, 1, N);

        theta  = [lambda; beta];
        thetaScaled = [lambdaScaled; betaScaled];
    elseif design==3
        % Covariance for the Gaussian copula
        corrCopula = [1 .4 .4; .4 1  0.4 ;  0.4 0.4 1];
        Z = mvnrnd([0 0 0], corrCopula, N); % draw from the copula
        U = normcdf(Z,0,1); % cdf
        % convert to data
        xMean = norminv(U(:,1),0,1)';
        lambda =   (unifinv(U(:,2),0.25, 0.75))';
        beta =  (1+(norminv(U(:,3),0,1)))';
        lambdaScaled = (lambda-.25)/sqrt(T)+0.5;
        betaScaled = (beta-1)/sqrt(T)+1;
        % Draw variance of the error term
        sigmaSq = exprnd(1, 1, N);
        % Draw data
        x = mvnrnd(xMean, eye(N), T); 
    end
    
    % Draw initial conditions
    y0 = randn([1, N]); % initial values, always independent
    % Draw idiosyncratic disturbances
    u = mvnrnd(zeros(N, 1), (diag(sigmaSq)), T);
    

end