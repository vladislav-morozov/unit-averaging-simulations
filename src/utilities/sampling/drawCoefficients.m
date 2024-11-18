% ===========================================================
% File: drawCoefficients.m
% Description: Draws a cross-sectional sample of coefficients and variances
% ===========================================================

function [thetaSample, sigmaSqSample, thetaComponentLabels] = ...
    drawCoefficients(approach, N, varNoise, seedCoefficients)
    % drawCoefficients Simulates coefficients and variance for a specified
    % model type.
    %
    % Generates AR(1) model coefficients `theta_{i1}` and `theta_{i2}` for
    % different data-generating process (DGP) approaches, along with error
    % variances `sigmaSqSample` for a given cross-sectional size.
    %
    % Args:
    %     approach (char): DGP for coefficients. Options are "unimodal", 
    %                      "bimodal", "bimodal_close", or "local".
    %     N (int): Cross-sectional sample size.
    %     varNoise (numeric): Scale for the variance of error term.
    %     seedCoefficients (int): Seed for random number generation.
    %
    % Returns:
    %     thetaSample (2xN matrix): Matrix of coefficients, where the first 
    %         row is `theta_{i1}` and the second row is `theta_{i2}` for 
    %         each cross-sectional unit.
    %     sigmaSqSample (1xN vector): Vector of variances for the error 
    %         term for each unit.
    %     thetaComponentLabels (Nx1 vector): Component label for each 
    %         observation, identifying mixture components (1 for unimodal).
    %
    % Example:
    %     [thetaSample, sigmaSqSample, thetaComponentLabels] = ...
    %         drawCoefficients("bimodal", 100, 0.5, 123);

    % Initialize random seed for reproducibility
    rng(seedCoefficients, 'philox');

    % Preallocate variables for coefficients
    lambdaSample = nan(1, N);
    betaSample = nan(1, N);
    % Default to all ones for unimodal distributions
    thetaComponentLabels = ones(N, 1);  

    % Select approach for coefficient generation
    switch approach
        case "unimodal"
            % "Unimodal": Draw AR(1) coefficients from a Beta distribution
            lambdaMean = 0.5;
            lambdaScale = 0.7;
            lambdaSample = lambdaMean + ...
                lambdaScale * (betarnd(5, 5, 1, N) - 0.5);

            % Draw theta_{i2} from normal distribution N(1, 1)
            betaSample = 1 + randn(1, N);

        case "bimodal"
            % "Bimodal": Draw coefficients from a mixture of two Beta
            % distributions
            lambdaCenters = [0.3, 0.65];
            lambdaScales = [0.4, 0.4];
            betaCenters = [0, 2];
            betaScales = [0.3, 0.3];
            for unitID = 1:N
                % Randomly select one of the two mixture components
                unitComponent = randi(2);  
                thetaComponentLabels(unitID) = unitComponent;
                lambdaSample(unitID) = lambdaCenters(unitComponent) + ...
                    lambdaScales(unitComponent) * (betarnd(5, 5) - 0.5);
                betaSample(unitID) = betaCenters(unitComponent) + ...
                    betaScales(unitComponent) * randn;
            end

        case "bimodal_close"
            % "Bimodal_close": Similar to "bimodal" but with closer Beta
            % mixture modes
            lambdaCenters = [0.39, 0.61];
            lambdaScales = [0.45, 0.45];
            betaCenters = [0, 2];
            betaScales = [0.3, 0.3];
            for unitID = 1:N
                % Randomly select one of the two mixture components
                unitComponent = randi(2);  
                thetaComponentLabels(unitID) = unitComponent;
                lambdaSample(unitID) = lambdaCenters(unitComponent) + ...
                    lambdaScales(unitComponent) * (betarnd(5, 5) - 0.5);
                betaSample(unitID) = betaCenters(unitComponent) + ...
                    betaScales(unitComponent) * randn;
            end

        case "local"
            % "Local": Localized parameters for AR(1) coefficients, using
            % the paper’s specifications for localized DGP.
            T = 60; % Time series length, fixed for "local" approach
            meanLambda = 0;
            meanBeta = 1;
            varianceBeta = 1;
            
            % Draw localized `theta_{i1}`
            etaLambda = 2 * (1 - abs(meanLambda)) * (2 * rand(1, N) - 1);
            lambdaSample = meanLambda + etaLambda / sqrt(T);
            
            % Draw localized `theta_{i2}`
            betaSample = meanBeta + ...
                sqrt(varianceBeta) * randn(1, N) / sqrt(T);

        otherwise
            error(['Invalid approach. Choose from "unimodal"', ...
                    '"bimodal", "bimodal_close", "local".'])
    end

    % Stack the coefficient vectors `lambdaSample` and `betaSample` into
    % `thetaSample`
    thetaSample = [lambdaSample; betaSample];

    % Draw variances for the error term from an exponential distribution
    sigmaSqSample = exprnd(varNoise, 1, N);
end
