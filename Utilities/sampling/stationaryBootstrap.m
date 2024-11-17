function [yBs, covarsBs] = stationaryBootstrap(y, covars, m, sampleLength)
% stationaryBootstrap Generates a bootstrap sample of a time series using
% the stationary bootstrap algorithm.
%
% This function implements the stationary bootstrap method proposed by
% Politis & Romano (1994). It generates a resampled version of the input
% time series data, preserving temporal dependencies through block
% sampling. Block lengths are determined probabilistically based on
% parameter m
%
% Args:
%     y (matrix): T x N matrix of time-series outcome variables. Rows
%         correspond to time points, and columns correspond to units.
%     covars (array): T x N x k array of covariates associated with the
%         time series. Second dimension indexes units, and third dimension 
%         indexes covariates.
%     m (scalar): Average block length parameter for the stationary
%         bootstrap algorithm. Determines the probability of starting a new
%         block (probability equal to 1/m).
%     sampleLength (scalar): Desired length of the output bootstrap sample.
%
% Returns:
%     yBs (matrix):  sampleLength x N  matrix of bootstrapped outcome 
%         variables. Rows correspond to sampled time points, and columns
%         correspond to units.
%     covarsBs (array): sampleLength x N x k array of bootstrapped 
%         covariates. Structure matches that of the input `covars`.
%
% Example:
%     >> y = [1; 2; 3; 4; 5; 6; 7; 8; 9; 10];
%     >> covars = rand(10, 1, 3);  % Example covariate data
%     >> [yBs, covarsBs] = stationaryBootstrap(y, covars, 4, 9);
%     >> disp(yBs);
%
% Reference:
%     Politis, D. N., & Romano, J. P. (1994). The Stationary Bootstrap. 
%     *Journal of the American Statistical Association*, 89(428), 1303–1313. 
%     DOI: 10.1080/01621459.1994.10476870.

    % Extract dimensions of the input covariates
    [T, N, k] = size(covars);

    % Probability of starting a new block
    accept = 1 / m;

    % Uniform random distribution for block decisions
    pd1 = makedist('Uniform');

    % Preallocate output arrays for bootstrap samples
    yBs = zeros(sampleLength, N);
    covarsBs = zeros(sampleLength, N, k);

    % Initialize sampling: draw a random starting index
    sampleIndex = randperm(T, 1);

    % Generate bootstrap sample
    for iSample = 1:sampleLength
        % Determine whether to continue the current block or start a new
        % one
        if random(pd1) >= accept
            % Continue the block: increment index
            sampleIndex = sampleIndex + 1;

            % Wrap around if the index exceeds the time series length
            if sampleIndex > T
                sampleIndex = 1;
            end
        else
            % Start a new block: draw a new starting index
            sampleIndex = randperm(T, 1);
        end

        % Add the selected data point to the bootstrap sample
        yBs(iSample, :) = y(sampleIndex, :);
        covarsBs(iSample, :, :) = covars(sampleIndex, :, :);
    end
end