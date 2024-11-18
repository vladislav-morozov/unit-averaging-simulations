function Q = ...
    uaObjectiveMatrix(estCoefs, estCovars, gradientEstimateTarget, ...
        targetIdx, unrestrictedBool)
    % uaObjectiveMatrix Constructs an objective matrix for optimal unit
    % averaging.
    %
    % This function builds a matrix used to minimize the mean-squared error
    % (MSE) in optimal unit averaging by constructing the Psi matrix that
    % encapsulates differences in coefficient estimates and associated
    % variances for unrestricted units.
    %
    % Args:
    %     estCoefs (matrix): k x N matrix of estimated coefficients, where 
    %         each column represents estimates from different 
    %         cross-sectional units.
    %     estCovars (array or cell array): Estimated variance-covariance 
    %         matrices, provided as either an N-element cell array or a 
    %         k x k x N array, where the third dimension indexes units.
    %     gradientEstimateTarget (vector): k-element gradient vector for 
    %         the target parameter with respect to estimated coefficients,
    %         evaluated at the target unit.
    %     targetIdx (int): Index of the target unit (1 ? targetIdx ? N).
    %     unrestrictedBool (vector, optional): Boolean vector indicating 
    %         unrestricted units. Defaults to all true (all units 
    %         unrestricted, fixed-N) if not provided.
    %
    % Returns:
    %     Q (matrix): Objective matrix representing the MSE minimization 
    %         problem, with the structure determined by the regime 
    %         (fixed-N or large-N).
    %
    % Example:
    %     Q = uaObjectiveMatrix(estCoefs, estCovars, ...
    %               gradientEstimateTarget, 1, [true, true, false]);

    % Default unrestrictedBool to all true if not provided
    if nargin < 5 || isempty(unrestrictedBool)
        unrestrictedBool = true(size(estCoefs, 2), 1);
    end

    % Number of unrestricted units
    numUnrestr = sum(unrestrictedBool);

    % Extract unrestricted coefficient estimates and covariances
    unrstrctCoefs = estCoefs(:, unrestrictedBool);
    unrstrctCovar = estCovars(unrestrictedBool);
    % Coefficient estimates for the target unit
    targetCoefEst = estCoefs(:, targetIdx);  

    % Number of coefficients per unit
    numCoef = size(estCoefs, 1); 

    % Compute differences in coefficient estimates from the target unit
    coefDif = unrstrctCoefs - targetCoefEst;
    coefDif = mat2cell(coefDif, numCoef, ones(1, numUnrestr));
    difMat = blkdiag(coefDif{:});

    % Construct the core Psi matrix using outer products of differences
    psiCore = difMat * ones(numUnrestr, numUnrestr) * difMat';

    % Add the diagonal variance-covariance matrices to Psi core
    if isa(unrstrctCovar, 'cell')
        % Unrestricted covariance provided as cell array
        varMat = unrstrctCovar;
    else
        % Unrestricted covariance provided as k x k x N array
        varMat = ...
            squeeze(mat2cell(unrstrctCovar, numCoef, numCoef, ...
                ones(1, numUnrestr)));
    end
    varMat = blkdiag(varMat{:});
    psiCore = psiCore + varMat;

    % Incorporate gradients to complete the Psi matrix
    gradientMatrix = kron(eye(numUnrestr), gradientEstimateTarget);
    Psi = gradientMatrix' * (psiCore + psiCore') * gradientMatrix;

    % Determine regime: large-N (restricted units present) or fixed-N (no
    % restricted units)
    if numUnrestr < size(estCoefs, 2)
        % Large-N regime: add elements for restricted units

        % Initialize Q matrix with an additional row and column for
        % restricted units
        Q = nan(numUnrestr + 1, numUnrestr + 1);
        b = nan(numUnrestr, 1);

        % Insert Psi matrix into Q
        Q(1:end-1, 1:end-1) = Psi;

        % Calculate and insert elements for the last row and column
        % (restricted units)
        mg = mean(estCoefs, 2);  % Mean of coefficient estimates
        for rowID = 1:numUnrestr
            % Calculate the cross-term component for restricted unit effects
            bRow = ...
                (unrstrctCoefs(:, rowID) - targetCoefEst) * ...
                (targetCoefEst - mg)';
            Q(rowID, end) = -gradientEstimateTarget' * bRow * ...
                gradientEstimateTarget;
            Q(end, rowID) = Q(rowID, end);  % Ensure symmetry
        end

        % Insert last element for the restricted units' interaction term
        Q(end, end) = (gradientEstimateTarget' * (targetCoefEst - mg))^2;

    else
        % Fixed-N regime: only unrestricted units present, so Q = Psi
        Q = Psi;
    end
end
