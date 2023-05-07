function [psiMatrix, psiOffDiag, psiDiag] = uaPsiMatrix(eta, V, d0, largeN, i0)
    % averagingPopulationPsi Construct the population psi matrix using
    % eta (k x N), the individual heterogeneity parameters, individual coefficient
    % covariance matrix, gradient d0. If largeN=1, it returns
    
    % If large-N options omitted, default to fixed-N regime
    if nargin<5
       largeN=0;
       i0=Inf;
    end
       
    if largeN==0 % Fixed-N regime
        [~, M] = size(eta);
        % Separately create bias elements that depend on eta (psiOffDiag)
        % and variance elemnts that go on the diagonal
        psiOffDiag = zeros(M, M);
        psiDiag = zeros(M, M);

        for i=1:M
            % Introduce bias terms
            for j=1:i
                psiOffDiag(i,j) = d0'*((eta(:, i)-eta(:,1))*(eta(:,j)-eta(:,1))'  )*d0 ;
            end
            % Add
            psiDiag(i, i) = d0'*( V(:, :,i)  )*d0;
        end
        psiOffDiag = psiOffDiag + triu(psiOffDiag',1); % add the above-diagonal part
        psiMatrix = psiOffDiag+ psiDiag; % combine
    else % Large-N regime
        % Follow the same strategy as above
        psiOffDiag = zeros(i0-1, i0-1);
        psiDiag = zeros(i0-1, i0-1);
        b = zeros(i0-1, 1);
        for i=1:i0-1
            % Introduce bias terms
            for j=1:i
                psiOffDiag(i,j) = d0'*((eta(:, i)-eta(:,1))*(eta(:,j)-eta(:,1))'  )*d0 ;
            end
            % Add
            psiDiag(i, i) = d0'*( V(:, :,i)  )*d0;
            b(i) = -d0'*(eta(:, i)-eta(:, 1))*(eta(:, 1))'*d0;
        end
        psiOffDiag = psiOffDiag + triu(psiOffDiag',1); % add the above-diagonal part
        psiMatrix = psiOffDiag+ psiDiag; % combine
        psiMatrix = [psiMatrix, b; b', (d0'*eta(:,1))^2];
        
    end

end
