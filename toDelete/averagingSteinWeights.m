function [weightsStein, psiStein] = averagingSteinWeights(eta, k, V, D)
    % averagingSteinWeights Finds optimal weights for combining the
    % individual and the MG estimators using large T fixed N MSEs

    [~, psiOffDiag, psiDiag] = averagingPopulationPsi(eta(:, 1:k), V, D);

    psiStein = [psiDiag(1, 1),  psiDiag(1, 1)/k; psiDiag(1, 1)/k,...
            sum(diag(psiDiag))/k^2 + sum(psiOffDiag, [1, 2])/k^2];
    
    options = optimoptions('quadprog','Display','none');
    weightsStein = quadprog(psiStein,zeros(2,1), [], [], ...
                ones(1, 2), 1,zeros(1, 2),[],[], options);
end