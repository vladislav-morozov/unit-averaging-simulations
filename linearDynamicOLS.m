function theta = linearDynamicOLS(y, x, T)
    % linearDynamicOLS Returns the OLS estimates for lambda_i and beta_i
    
    theta = [y(1:T), x]\y(2:T+1);
end