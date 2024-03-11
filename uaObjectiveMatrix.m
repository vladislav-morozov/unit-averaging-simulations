function Q = uaObjectiveMatrix(est_coefs, est_covars, gradient_estimate_target, target_idx, unrestricted_bool)
% Builds the objective matrix for unit averaging. Accommodates both
% finite-N and large-N approximations.
% IMPORTANT: in contrast to the paper, the estimators and the covariances are
% not multiplied by T. One should use directly use the covariances returned
% by the estimation algorithm
%
% Args:
%     unrestricted_bool: None, boolean array or vector with dimension
%                        matching est_coefs. True means the unit is
%                        unrestricted. If None, all units are unrestricted

if nargin < 5 || isempty(unrestricted_bool)
    unrestricted_bool = true(size(est_coefs));
end

% Compute the matrix on unrestricted units
unrstrct_coefs = est_coefs(unrestricted_bool);
unrstrct_covar = est_covars(unrestricted_bool);

Psi = zeros(length(unrstrct_coefs), length(unrstrct_coefs));

for i = 1:length(unrstrct_coefs)
    for j = 1:length(unrstrct_coefs)
        
        % Difference between estimates
        coef_dif_i = unrstrct_coefs(i) - est_coefs(target_idx);
        coef_dif_j = unrstrct_coefs(j) - est_coefs(target_idx);
        psi_ij = coef_dif_i * coef_dif_j';
        % add covariance when appropriate
        if i == j
            psi_ij = psi_ij + unrstrct_covar(i);
        end
        % Multiply by gradient
        psi_ij = gradient_estimate_target * psi_ij * gradient_estimate_target';
        
        % Set the corresponding element
        Psi(i, j) = psi_ij;
    end
end

% If in large-N regime, add the outer row and column of restricted units
if sum(unrestricted_bool) < length(est_coefs)
    Q = zeros(length(unrstrct_coefs) + 1, length(unrstrct_coefs) + 1);
    Q(1:end-1, 1:end-1) = Psi;
    b = zeros(length(unrstrct_coefs), 1);
    
    % Fill out the elements of the b vector
    mg = mean(est_coefs);
    for i = 1:length(b)
        b_i = (unrstrct_coefs(i) - est_coefs(target_idx)) * (est_coefs(target_idx) - mg)';
        Q(i, end) = -gradient_estimate_target * b_i * gradient_estimate_target';
        Q(end, i) = Q(i, end);
    end
    
    % Insert the last element
    Q(end, end) = (gradient_estimate_target * (est_coefs(target_idx) - mg))^2;
    
else
    Q = Psi;
end

% Return the Q matrix
end
