function Q = ...
    uaObjectiveMatrix(estCoefs, estCovars, gradientEstimateTarget, targetIdx, unrestrictedBool)
% uaObjectiveMatrix Builds the objective matrix for unit averaging. 
% Accommodates both finite-N and large-N approximations with the
% unrestrictedBool argument
%
% Args:
%     estCoefs -- kxN matrix, columns index different units
%     estCovars -- N-cell array where elements are covariance matrices or 
%                  kxkxN matrix, third index indices units
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
unrstrctCovar = estCovars(unrestrictedBool);
targetCoefEst = estCoefs(:, targetIdx);

% Allocate Psi matrix
Psi = nan(numUnrestr, numUnrestr);

% % Fill out the Psi matrix
% for colID = 1:numUnrestr
%     for rowID = 1:colID
%     
%         
%         % Difference between estimates
%         coefDifRow = unrstrctCoefs(:, rowID) - targetCoefEst;
%         coefDifCol = unrstrctCoefs(:, colID) - targetCoefEst;
%         psiRowCol = coefDifRow * coefDifCol';
%         
%         % add covariance when appropriate
%         if rowID == colID
%             psiRowCol = psiRowCol + squeeze(unrstrctCovar(:, :, rowID));
%         end
%         % Multiply by gradient
%         psiRowCol = ...
%             gradientEstimateTarget'*psiRowCol*gradientEstimateTarget;
%         
%         % Set the corresponding element
%         Psi(rowID, colID) = psiRowCol;
%         Psi(colID, rowID) = psiRowCol;
%     end
% end



% Experimental
numCoef = size(estCoefs, 1);
% psiCore = zer?os(numUnrestr*numCoef, numUnrestr*numCoef);

% Evaluate the outer products of differences
coefDif = unrstrctCoefs - targetCoefEst;
coefDif = mat2cell(coefDif, numCoef, ones(1, numUnrestr));
difMat = blkdiag(coefDif{:});
psiCore = difMat*ones(numUnrestr,numUnrestr)*difMat';

% Add a diagonal of covariances
if class(unrstrctCovar) ~= "cell"
    varMat = squeeze(mat2cell(unrstrctCovar, numCoef,numCoef, ones(1,numUnrestr)));
else
    varMat = unrstrctCovar;
end
varMat = blkdiag(varMat{:});
psiCore = psiCore+varMat;

% for colID = 1:numUnrestr
%     for rowID = 1:colID
%     
%         
%         Difference between estimates
%         coefDifRow = unrstrctCoefs(:, rowID) - targetCoefEst;
%         coefDifCol = unrstrctCoefs(:, colID) - targetCoefEst;
%         psiRowCol = coefDifRow * coefDifCol';
%         
%         add covariance when appropriate
%         if rowID == colID
%             psiRowCol = ...
%                 0.5*(psiRowCol + squeeze(unrstrctCovar(:, :, rowID)));
%         end
%         
%         
%         Set the corresponding element
%         posRow = (rowID-1)*numCoef+1;
%         posCol = (colID-1)*numCoef+1;
%         psiCore(posRow:posRow+numCoef-1, posCol:posCol+numCoef-1) = psiRowCol;
%         psiCore(posCol:posCol+numCoef-1, posRow:posRow+numCoef-1) = psiRowCol';
%     end
% end 
gradientMatrix = kron(eye(numUnrestr), gradientEstimateTarget);
Psi = gradientMatrix'*(psiCore+psiCore')*gradientMatrix;


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
