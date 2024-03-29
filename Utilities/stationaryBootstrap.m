function [yBs, covarsBs] = stationaryBootstrap(y, covars, m, sampleLength)
% STATIONARY BOOTSTRAP Returns a bootstraped sample of the time-series
% "data" of length "sampleLength. The algorithm used is stationary bootstrap
% from 1994 Politis & Romano.
%    sample = StationaryBootstrap(data, m, sampleLength)
%
% Input:
%    data = n x 1 vector containing the time-series that will be
%                       bootstrapped
%    m = 1 x 1 parameter m of the stationary bootstrap algorithm indicating
%        the average length of blocks in the final sample
%    sampleLength = 1 x 1 integer setting the lenght of the output sample
%
% Output:
%    sample =       sampleLength x 1 vector containing the bootstraped sample
%
% Example:
%    >> data = [1;2;3;4;5;6;7;8;9;10]
%    >> StationaryBootstrap(data,4,9)
%    ans = [9; 10; 1; 2; 5; 6; 5; 6; 7]
%
% Original source:
%     Dimitris N. Politis & Joseph P. Romano (1994) The Stationary Bootstrap, Journal of the American Statistical
%     Association, 89:428, 1303-1313, DOI: 10.1080/01621459.1994.10476870
%

% Extract data dimensions
[T, N, k] = size(covars);

% Prepare for sampling
accept = 1/m;
pd1 = makedist('Uniform'); 

% Allocate output arrays
yBs = zeros(sampleLength, N);
covarsBs = zeros(sampleLength, N, k);

% Initialize the sampler
sampleIndex = randperm(T, 1);
for iSample = 1:sampleLength
    % Check whether the block stops or not
    if random(pd1)>=accept
        % Block continues because step is accepted
        
        % Increment index
        sampleIndex = sampleIndex+1;
        % If the block arrives at the end of the data, loop the idx around
        % to the start of the date
        if sampleIndex > T
            sampleIndex = 1;
        end
    else
        % Step is rejected, therefore a new block is started
        
        % Draw new starting pointx
        sampleIndex = randperm(T, 1);
    end
    % Add the selected datapoint to the sample
    yBs(iSample, :) = y(sampleIndex, :); 
    covarsBs(iSample, :, :) = covars(sampleIndex, :, :);
end
end
