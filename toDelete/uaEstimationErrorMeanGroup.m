function errEstMG =...
            uaEstimationErrorMeanGroup(Nvalues, targetValue,unitEst)
% uaEstimationErrorMeanGroup Evaluates performance of unit
% averaging with mean group weights
% 
% Arguments: 
% 1. Nvalues -- vector of values for N to consider
% 2. targetValue -- target value to estimate
% 3. unitEst -- individual estimators, units arranged as columns
%
% Returns   vector of errors, size of Nvalues
% 1. errEstMG -- errors with MG weights
 
    Nlen = length(Nvalues);
    errEstMG = zeros(1, Nlen);
 
    
    for i=1:Nlen
        Ncurrent = Nvalues(i); % select N
        avgEst  = mean(unitEst(1:Ncurrent));  % predict
        errEstMG(i) = avgEst - targetValue;  % tabulate error
    end

end