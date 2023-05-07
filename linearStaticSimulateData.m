function [y, x, u] = linearStaticSimulateData(beta, sigmaSq, T, design, seedData)
    % linearDynamicSimulateData Simulates data for the first MC experiment
    % The model is y_{it} = \lambda_{i} y_{it-1} + \beta_{i}x_{it} + u_{it}
    % Draws data according to design chosen and returns data and parameter
    % vectors. Designs values accepted -- '1', '2'
    
    [~, N] = size(beta);
    
    if design==1 % Draw data
        rng(seedData,'philox')
        x = mvnrnd(ones(1, N), eye(N), T); 
        x(:, :, 2) = mvnrnd(ones(1, N), eye(N), T); 
        u = mvnrnd(zeros(N, 1), (diag(sigmaSq)), T);

    end
    
    
    % Simulate data
    y = x(:,:,1).*repmat(beta(1,:), T,1) + x(:,:,2).*repmat(beta(2,:), T,1)+u; 
    
end