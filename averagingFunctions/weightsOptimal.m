function optWeights = ...
    weightsOptimal(estCoefs, estCovars, gradientEstimateTarget, ...
        targetID, unrestrictedBool)
    % weightsOptimal Computes optimal unit averaging weights based on
    % estimated coefficients and covariances.
    %
    % This function calculates the optimal weights for averaging across
    % units to minimize mean-squared error (MSE), supporting both fixed-N
    % and large-N approaches for unrestricted units.
    %
    % Args:
    %     estCoefs (matrix): k x N matrix of estimated coefficients, where 
    %         each column represents estimates from different 
    %         cross-sectional units.
    %     estCovars (3D array): k x k x N array of estimated 
    %         variance-covariance matrices, where the third dimension
    %         indexes cross-sectional units.
    %     gradientEstimateTarget (vector): k-element gradient vector for 
    %         the target parameter with respect to estimated coefficients,
    %         evaluated at the target unit.
    %     targetID (int): Index of the target unit (1 ? targetID ? N).
    %     unrestrictedBool (vector, optional): Boolean vector indicating 
    %         unrestricted units. Defaults to all true 
    %         (all units unrestricted, fixed-N) if not provided.
    %
    % Returns:
    %     optWeights (vector): N-element vector of optimal averaging
    %          weights.
    %
    % Example:
    %     optWeights = weightsOptimal(estCoefs, estCovars, ...
    %       gradientEstimateTarget, 1, [true, true, false]);

    % Default unrestrictedBool to all true if not provided (fixed-N regime)
    if nargin < 5 || isempty(unrestrictedBool)
        unrestrictedBool = true(size(estCoefs, 2), 1);
    end

    % Construct the objective matrix using uaObjectiveMatrix function
    Q = ...
        uaObjectiveMatrix(estCoefs, estCovars, gradientEstimateTarget,...
               targetID, unrestrictedBool);
    numCoords = size(Q, 1);  % Determine size of Q

    % Symmetrize Q to prevent numerical warnings in quadprog
    Q = (Q + Q') / 2;

    % Set up optimization options for the quadratic programming solver
    options = optimoptions('quadprog', 'Display', 'none');

    % Solve the quadratic programming problem to obtain optimal weights
    % Subject to: sum(weights) = 1 and weights >= 0
    xOpt = quadprog(Q, zeros(numCoords, 1), [], [], ones(1, numCoords), 1, ...
                    zeros(1, numCoords), [], [], options);

    % Allocate the output vector for weights across all units
    optWeights = nan(size(estCoefs, 2), 1);

    % Determine the number of restricted units
    numRestr = size(estCoefs, 2) - sum(unrestrictedBool);

    if numRestr == 0
        % Fixed-N regime: all units unrestricted, directly assign xOpt to
        % optWeights
        optWeights = xOpt;
    else
        % Large-N regime: some units are restricted Allocate weights for
        % unrestricted units and distribute remaining weight across
        % restricted units

        % Insert weights for unrestricted units
        unresWeights = xOpt(1:end-1);
        optWeights(unrestrictedBool) = unresWeights;

        % Calculate and assign the split weight for restricted units
        resWeightInd = xOpt(end) / numRestr;
        optWeights(~unrestrictedBool) = resWeightInd;
    end
end
