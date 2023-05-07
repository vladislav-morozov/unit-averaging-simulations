%% Parameter block
% Define parameters used in all simulations here

% Data dimensions
N= 26;
T = 100;


% For Dynamic model: select default mean and variance parameters for
% coefficients
meanBeta = 1; % mean for beta coefficient
varianceBeta = 1; % variance for beta
meanCoef = [0; meanBeta]; % define the mean vector for coefficients

% Number of replications
numReplications = 10000;

% Negative weights, set to 0 to enforce nonnegativity, 1 to allow
nW = 0;

% Simulation design choice, options are 1 and 2
design =1;

% Select a seed to drawing coefficients
seedCoefficients = 15;


% Draw coefficients and variances
betaOriginal = ones(2, N);
        
% Set range to be explored
betaOriginal(1,:) = (-12.5:1:12.5)/26;


etaTrue = betaOriginal - repmat(meanCoef, 1, N); % obtain deviations from the mean
betaScaled = meanCoef+etaTrue/sqrt(T);

% Create objects for later use
 % Create equal weights
 equalWeightsAll = ones(N,1)/(N);



%% Evolution of weights
 
 
 weightsSample2 = zeros(N, N);
 for j=1:N
  
     
   % Choose gradient, gradient describes which parameter is estimated
    D = [1; 0];% lambda_1
%    D = [betaScaled(2, 1)/(1-betaScaled(1, 1))^2; 1/(1-betaScaled(1, 1))]; % long-run
   % D = [2*betaScaled(2,1); 2*betaScaled(1, 1)] % norm
% 
%    Compute variance and population Psi matrices
%    V = linearDynamicTrueCovariance(betaScaled(1,:), betaScaled(2,:),sigmaSq, N,design);    % compute population-fixed variances
%    [psiMatrixTrue, ~, ~] = averagingPopulationPsiPosition(etaTrue, V, D, j); % recompute population Psi
%    Special computations for Stein-type estimator 
%  
%    
   % Recreate temporary variance vectors
   weightsSampleInLoop = zeros(N, numReplications);

   % Inner loop: drawing data
    parfor i=1:numReplications
        [~,~ , sigmaSq] = linearDynamicDrawCoefficients(N, T, design,... 
            meanBeta, varianceBeta,   i);
        % Draw data, estimate coefficients and the gradent
        [y, x, u] = linearDynamicSimulateData(betaScaled, sigmaSq, T, design, i); % draw data
        betaHat = linearStaticEstimators(y, x);% estimate coefficients
        etaEstimated = sqrt(T)*(betaHat - repmat(mean(betaHat,2),1, N));
        Vest = linearDynamicVarianceEstimator(y, x, betaHat); % estimate variance
        Dest = D;
        
        % Compute Psi matrics based on data
        psiMatrixSample = averagingSamplePsiPosition(betaHat*sqrt(T), Vest, Dest, j);% estimate sample psi

        weightsSampleInLoop(:,i)= averagingOptimalWeights(psiMatrixSample(1:N,1:N), nW);
        [j, i]
    end
    
    weightsSample2(:, j) = mean(weightsSampleInLoop,2);
    
 end
 

 %%
 
 % Weight distribution
 figure
 area(betaOriginal(1,:), weightsSample2')

xlabel('\eta_1')
title('Proportions of units in averaging estimators as unit of interest changes, estimating \lambda_i, N=25')
xlim([min(betaOriginal(1,:)), max(betaOriginal(1,:)) ] )
ylim([0, 1])
xlabel('i')
ylabel('Weights')

% Evolution of weights
figure
yLabels = {'Unit 1', 'Unit 6', 'Unit 11', 'Unit 16', 'Unit 21', 'Unit 26'};
stackedplot(betaOriginal(1,:), weightsSample2([1,6,11,16, 21, 26], :)','DisplayLabels',yLabels, 'Marker','o')
xlabel('i')

title('Weights of selected units as unit of interest changes, estimating \lambda_i, N=25')

% Own weight
figure
plot(betaOriginal(1,:), diag(weightsSample2), 'Marker','o');
xlabel('i')
ylabel('w_i')
title('Own weight of unit i, estimating \lambda_i, N=25')