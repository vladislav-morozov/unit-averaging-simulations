function Q = ...
    uaObjectiveMatrix(estCoefs, estCovars, gradientEstimateTarget, targetIdx, unrestrictedBool)
% uaObjectiveMatrix Builds the objective matrix for unit averaging. 
% Accommodates both finite-N and large-N approximations with the
% unrestrictedBool argument
%
% Args:
%     estCoefs -- kxN matrix, columns index different units
%     estCovars -- kxkxN matrix, third index indices units
%     targetIdx -- positive integer <= N, id of the target unit
%     gradientEstimateTarget -- k-vector; gradient of the target parameter
%                               with respect to thetas, evaluate at the
%                               target unit
%     unrestricted_bool: (optional) empty or Boolean vector with dimension
%                        matching estCoefs. True means the unit is
%                        unrestricted. If None, all units are unrestricted
%
% Returns:
%     Q -- matrix describing the MSE minimization problem

 
% If no vector of restricted units is supplied, assume fixed-N behavior
if nargin < 5 || isempty(unrestrictedBool)
    unrestrictedBool = true(size(estCoefs, 2), 1);
end

% Obtain number of unrestricted units
numUnrestr = sum(unrestrictedBool);

% Extract unrestricted coefficients and covariates
unrstrctCoefs = estCoefs(:, unrestrictedBool);
unrstrctCovar = estCovars(:, :, unrestrictedBool);
targetCoefEst = estCoefs(:, targetIdx);

% Allocate Psi matrix
Psi = nan(numUnrestr, numUnrestr);

% Fill out the Psi matrix
for rowID = 1:numUnrestr
    for colID = 1:numUnrestr
        
        % Difference between estimates
        coefDifRow = unrstrctCoefs(:, rowID) - targetCoefEst;
        coefDifCol = unrstrctCoefs(:, colID) - targetCoefEst;
        psiRowCol = coefDifRow * coefDifCol';
        
        % add covariance when appropriate
        if rowID == colID
            psiRowCol = psiRowCol + squeeze(unrstrctCovar(:, :, rowID));
        end
        % Multiply by gradient
        psiRowCol = ...
            gradientEstimateTarget'*psiRowCol*gradientEstimateTarget;
        
        % Set the corresponding element
        Psi(rowID, colID) = psiRowCol;
    end
end

% If in large-N regime, add the outer row and column of restricted units
if numUnrestr < size(estCoefs, 2)
    % Allocate Q and b
    Q = nan(numUnrestr + 1, numUnrestr + 1);
    b = nan(numUnrestr, 1);
    
    % Insert Psi ino Q
    Q(1:end-1, 1:end-1) = Psi;
    
    % Fill out the elements of the b vector
    mg = mean(estCoefs, 2);
    for rowID = 1:length(b)
        bRow = (unrstrctCoefs(:, rowID) - targetCoefEst) * (targetCoefEst - mg)';
        Q(rowID, end) = -gradientEstimateTarget' * bRow * gradientEstimateTarget;
        Q(end, rowID) = Q(rowID, end);
    end
    
    % Insert the last element
    Q(end, end) = (gradientEstimateTarget' * (targetCoefEst - mg))^2;
    
else
    Q = Psi;
end

% Return the Q matrix
end
