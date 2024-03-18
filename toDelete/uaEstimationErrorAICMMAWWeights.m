function [errEstAIC, errEstMMA] =...
            uaEstimationErrorAICMMAWWeights(thetaHat, y, x, Nvalues, targetValue, unitEst)
% uaEstimationErrorAICMMAWWeights Evaluates performance of unit
% averaging with exponential AIC and MMA weights weights.
% 
% Arguments: 
% 1. Parameter estimates, units arranged as columns
% 2, 3. Data y, x
% 4. Nvalues -- vector of values for N to consider
% 5. targetValue -- target value to estimate
% 6. unitEst -- individual estimators, units arranged as columns
%
% Returns   vectors of errors, size of Nvalues
% 1. errEstAIC -- errors with exponential AIC weights
% 2. errEstMMA -- errors with MMA weights (observe that for unit averaging
%                 with the same unit models this is equivalent to mean AIC/BIC model choice)
    
    Nlen = length(Nvalues);
    errEstAIC = zeros(1, Nlen);
    errEstMMA = zeros(1, Nlen);
    
    for i=1:Nlen
        Ncurrent = Nvalues(i); % select N
        [aicW, mmaW] = uaLinearAICMMweights(thetaHat, y, x, Ncurrent);
         
        % AIC
        avgEst  = unitEst(1:Ncurrent)*aicW;  % predict
        errEstAIC(i) = avgEst - targetValue;  % tabulate error
        
        % MMA  % predict
        avgEst  = unitEst(1:Ncurrent)*mmaW;
        errEstMMA(i) = avgEst - targetValue;  % tabulate error% N
        
    end

end