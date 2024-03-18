function psiMatrix = averagingSamplePsiPosition(individualEstimators,varianceEstimators,D, pos)

    [~, N] = size(individualEstimators);
    psiMatrix = zeros(N, N);
    eta = individualEstimators;
%     thetaMG = mean(eta, 2);
    V = varianceEstimators;
%     psiMatrix(1,1) = D'*(eta(:,1)-thetaMG)*(eta(:,1)-thetaMG)'*D;
    for i=1:N
        for j=1:N
        psiMatrix(i,j) = D'*((eta(:, i)-eta(:,pos))*(eta(:,j)-eta(:,pos))'+ V(:, :,i)*(i==j)  )*D ; 
        end
    end
%     psiMatrix = psiMatrix + triu(psiMatrix', 1);
end