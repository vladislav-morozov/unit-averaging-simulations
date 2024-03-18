function [errEstFixed, errEstLarge, errEstLargeFavorable, errEstLargeUnfavorable] =...
            uaEstimationErrorPluginWeights(etaEst, V, d1, i0, Nvalues, trueValue, unitEst)
% uaEstimationErrorInfeasibleWeights Evaluates performance of unit
% averaging with infeasible weights. Arguments: 
% 1. Estimated heterogeneity parameters etaEst, units arranged as columns
% 2. Population asymptotic variances V, units arranged along the third
%    dimension
% 3. d1 -- gradient of parameter of interest
% 4. i0 -- cutoff parameter for non-negligible weights
% 5. Nvalues -- vector of values for N to consider
% 6. trueValue -- target value to estimate
% 7. unitEst -- individual estimators, units arranged as columns
% Returns four vectors of errors, size of Nvalues
% 1. errEstFixed -- applying fixed-N regimex
% 2. errEstLarge -- applying large-N regime in the order supplied
% 3. errEstLargeFavorable -- large-N regime, units ordered by least bias
% 4. errEstLargeUnfavorable -- large-N reegime, units ordered by most bias
    
    Nlen = length(Nvalues);
    errEstFixed = zeros(1, Nlen);
    errEstLarge = zeros(1, Nlen);
    errEstLargeFavorable = zeros(1, Nlen);
    errEstLargeUnfavorable = zeros(1, Nlen);
    
    psiFixed = uaPsiMatrix(etaEst, V, d1);    % for fixed N mode the Psi matrix can be computed once
    for i=1:Nlen
        Ncurrent = Nvalues(i); % select N
        % fixed-N approximations
        psiLoop = psiFixed(1:Ncurrent, 1:Ncurrent);    % slice out Psi
        weightLoop = uaOptimalWeights(psiLoop);     % compute weights
        avgEst  = unitEst(1:Ncurrent)*weightLoop;  % predict
        errEstFixed(i) = avgEst - trueValue;  % tabulate error
        
        % large-N approximations: random. Since units are not ordered, we
        % can just pick the first ones for this case
        if Ncurrent<i0 % switch to fixed-N mode in case i0 is larger than current sample
            errEstLarge(i) = errEstFixed(i);    
        else % i_0<Ncurrent
            psiLoop =uaPsiMatrix(etaEst, V, d1, 1, i0); % it's convenient to construct the Q matrix for every value of
            weightLoop = uaOptimalWeights(psiLoop);     % compute weights
            avgEst  = uaAvgEstLargeN(unitEst, weightLoop);  % predict
            errEstLarge(i) = avgEst - trueValue;  % tabulate error% N
        end

%         large-N: favorable bias
        if Ncurrent<i0 % switch to fixed-N mode in case i0 is larger than current sample
            errEstLargeFavorable(i) = errEstFixed(i);    
        else % i_0<Ncurrent
            % Obtain units with least bias
            unitsSel = unitEst(1:Ncurrent);
            dists = abs(unitsSel-trueValue);
            [~, I] = sort(dists);
            psiLoop =uaPsiMatrix(etaEst(:, I), V(:, :, I), d1, 1, i0); % it's convenient to construct the Q matrix for every value of
            weightLoop = uaOptimalWeights(psiLoop);     % compute weights
            avgEst  = uaAvgEstLargeN(unitsSel(I), weightLoop);  % predict
            errEstLargeFavorable(i) = avgEst - trueValue;  % tabulate error% N
        end
%         large-N: least favorable bias
        if Ncurrent<i0 % switch to fixed-N mode in case i0 is larger than current sample
            errEstLargeUnfavorable(i) = errEstFixed(i);    
        else % i_0<Ncurrent
            % Obtain units with least bias
            unitsSel = unitEst(1:Ncurrent);
            dists = abs(unitsSel-trueValue);
            [~, I] = sort(dists, 'descend');
            I = circshift(I, 1); % ensure that 1 is always in the first position
            psiLoop =uaPsiMatrix(etaEst(:, I), V(:, :, I), d1, 1, i0); % it's convenient to construct the Q matrix for every value of
            weightLoop = uaOptimalWeights(psiLoop);     % compute weights
            avgEst  = uaAvgEstLargeN(unitsSel(I), weightLoop);  % predict
            errEstLargeUnfavorable(i) = avgEst - trueValue;  % tabulate error% N
        end
    end


end