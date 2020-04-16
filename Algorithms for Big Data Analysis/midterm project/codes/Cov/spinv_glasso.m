function [X, out] = spinv_glasso(S, rho, maxIter, epsilon)
    n = size(S,1);
    if nargin < 4
        epsilon = 1e-12; 
    end
    if nargin < 3
        maxIter = 100; 
    end
    W = S + rho * eye(n);
    W0 = W; beta = zeros(n-1,n);
    shr = sum(sum(abs(S)))-sum(abs(diag(S)));
    i = 0;
    % Graphical Lasso loop
    while i < maxIter
        i = i+1;
        for j = n:-1:1
            ind = setdiff(1:n,j);
            opts = [];
            s = l1_fprox_nesterov(zeros(n-1,1),W(ind,ind), S(ind,j), rho, opts);%use proximal gradient method
            beta(:,j) = s;
            W(ind,j) = W(ind,ind) * beta(:,j);
            W(j,ind) = W(ind,j)';
        end
        if sum(sum(abs(W-W0))) < epsilon*shr
            break;
        end
        W0 = W;
    end
    if i == maxIter
        fprintf('Maximum number of iteration reached.\n');
    end
    X = zeros(n);
    for j = 1 : n
        ind = setdiff(1 : n, j);
        X(j,j) = 1 / (W(j,j) - W(ind,j)' * beta(:,j));
        X(ind, j) = - X(j, j) * beta(:,j);
    end
    X = (X + X')/2;

    out = [];
    out.dualgap =  - n + trace(S*X)  + rho*norm(vec(X), 1);
    out.iter = i;
    out.value = log(det(X)) - trace(S*X) - rho * norm(vec(X), 1);
end