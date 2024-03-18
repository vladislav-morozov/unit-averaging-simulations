function [psiMatrix, psiOffDiag, psiDiag] = averagingPopulationPsi(eta, V, D)
    % averagingPopulationPsi Construct the population psi matrix
    
        
    [~, N] = size(eta);
    psiOffDiag = zeros(N, N);
    psiDiag = zeros(N, N);
%     psiMatrix(1,1) = D'*eta(:,1)*eta(:,1)'*D;
    for i=1:N
        for j=1:i
        psiOffDiag(i,j) = D'*((eta(:, i)-eta(:,1))*(eta(:,j)-eta(:,1))'  )*D ; 
        end
        psiDiag(i, i) = D'*( V(:, :,i)  )*D;
    end
    psiOffDiag = psiOffDiag + triu(psiOffDiag',1);
    psiMatrix = psiOffDiag+ psiDiag;
end
