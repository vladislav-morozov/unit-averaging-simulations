function opt_weights = uaAveragingWeights(est_coefs, est_covars, gradient_estimate_target, target_idx, unrestricted_bool)
% Checking for fixed-N/large-N regime
% If no value is supplied, fixed-N regime is the default
if nargin < 5 || isempty(unrestricted_bool)
    unrestricted_bool = true(size(est_coefs));
end

if length(unrestricted_bool) ~= length(est_coefs)
    error('ValueError')
end

% Check for missing data
non_missing_bool = true(size(est_coefs));
for ind_coef_id = 1:length(est_coefs)
    ind_coef = est_coefs(ind_coef_id);
    if isempty(ind_coef)
        non_missing_bool(ind_coef_id) = false;
    end
end

% Adjust the target index
target_idx = sum(~non_missing_bool(1:target_idx));

% Construct the objective function
Q = ua_objective_matrix(est_coefs(non_missing_bool), est_covars(non_missing_bool), gradient_estimate_target, target_idx, unrestricted_bool(non_missing_bool));
num_coords = size(Q, 1);
Q = matrix(Q);
lin_term = matrix(zeros(num_coords, 1));

% Specify the constraints
ineq_mat = matrix(-eye(num_coords));
ineq_vec = matrix(zeros(num_coords, 1));
eq_lhs = matrix(ones(1, num_coords));
eq_rhs = matrix(1.0);

% Minimize
options.show_progress = false;
q_res = qp(Q, lin_term, ineq_mat, ineq_vec, eq_lhs, eq_rhs);

% Fill the weights to the non-missing units, missing units receive zero
opt_weights = zeros(size(est_coefs));
x_opt = squeeze(q_res.x);

num_restr = length(est_coefs) - sum(unrestricted_bool);
if num_restr == 0
    % In the fixed-N regime sufficient to just spread the non-missing weights
    opt_weights(non_missing_bool) = x_opt;
else
    % If in large-N mode, we need to correctly allocate the weights
    % into correct non-missing and non-restricted positions first, and then
    % split the residual weight across restricted units
    % Insert unrestricted weights
    unres_weights = x_opt(1:end-1);
    opt_weights(non_missing_bool & unrestricted_bool) = unres_weights;
    % Split restricted weights
    res_weight_ind = x_opt(end) / num_restr;
    opt_weights(~unrestricted_bool) = res_weight_ind;
end

% Return the opt_weights
end
