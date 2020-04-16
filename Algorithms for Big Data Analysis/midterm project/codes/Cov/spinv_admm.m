function [X, out] = spinv_admm(S, sigma, opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% modified from the paper
%  An efficient admm algorithm for high dimensional precisionmatrix estimation via penalized quadratic loss.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 3
    opts = [];
end
% innitial mu
if ~isfield(opts,'mu_0'); opts.mu_0 = 0.1; end
% stopping criterion
if ~isfield(opts,'tol_rel'); opts.tol_rel = 1e-12; end 
if ~isfield(opts,'tol_gap'); opts.tol_gap = 1e-6; end
if ~isfield(opts,'maxIter'); opts.maxIter = 5e3; end
n = size(S,1);
% initial X: sparsest, initial Y: densest
X = zeros(n,n); Y = ones(n,n); Lambda = zeros(n,n);
mu = opts.mu_0;
dualgap = inf;
S2 = S*S; % store S^2
f1 = inf; iter = 0; maxIter = opts.maxIter;
while iter < maxIter
    iter = iter + 1;
    % update X
    Xp = X;
    B = mu*sigma/2*S2; A = B + eye(n);
    C = mu*sigma*S + Y - Lambda;
    X = sylvester(A,B,C);
    % update Y
    Yp = Y;
    Y0 = X + Lambda;
    Y = sign(Y0).*max(abs(Y0)-mu,0);
    Lambda = X - Y + Lambda;
    % check stopping
    SX = S*X; SY = S*Y;
    dualgap_X = norm(vec(X),1) + sigma*(norm(SX,'fro')^2-trace(SX));
    dualgap_Y = norm(vec(Y),1) + sigma*(norm(SY,'fro')^2-trace(SY));
    dualgap = max(abs([dualgap_X, dualgap_Y]));
    % update value
    f0 = f1;
    f1 = sigma*norm(SX,'fro')^2/2 + norm(vec(Y),1);

    frel = abs(f1-f0)/max(abs([f1,f0,1]));
    Xrel = norm(X-Xp,'fro')/max([1,norm(X,'fro'),norm(Xp,'fro')]);
    Yrel = norm(Y-Yp,'fro')/max([1,norm(Y,'fro'),norm(Yp,'fro')]);
    
    if (max([frel, Xrel, Yrel]) < opts.tol_rel) || dualgap < opts.tol_gap
        break
    end
end

try
    chol(X+eps*eye(size(X)));
catch %continue iteration with projection
    while iter < maxIter
    iter = iter + 1;
    % update X
    Xp = X;
    B = mu*sigma/2*S2; A = B + eye(n);
    C = mu*sigma*S + Y - Lambda;
    X = sylvester(A,B,C);
    try
        chol(X+eps*eye(size(X)));
    catch
        [V,D] = eig(X);
        d = diag(D); d = max(d,0);
        X = (V*diag(d))*V';
    end
    % update Y
    Yp = Y;
    Y0 = X + Lambda;
    Y = sign(Y0).*max(abs(Y0)-mu,0);
    Lambda = X - Y + Lambda;
    % check stopping
    SX = S*X; SY = S*Y;
    dualgap_X = norm(vec(X),1) + sigma*(norm(SX,'fro')^2-trace(SX));
    dualgap_Y = norm(vec(Y),1) + sigma*(norm(SY,'fro')^2-trace(SY));
    dualgap = max(abs([dualgap_X, dualgap_Y]));
    % update value
    f0 = f1;
    f1 = sigma*norm(SX,'fro')^2/2 + norm(vec(Y),1);

    frel = abs(f1-f0)/max(abs([f1,f0,1]));
    Xrel = norm(X-Xp,'fro')/max([1,norm(X,'fro'),norm(Xp,'fro')]);
    Yrel = norm(Y-Yp,'fro')/max([1,norm(Y,'fro'),norm(Yp,'fro')]);

    if (max([frel, Xrel, Yrel]) < opts.tol_rel) || dualgap < opts.tol_gap
        break
    end
    end
end

out = [];
out.iter = iter; 
if iter == opts.maxIter
    fprintf('maximum number of iteration achieved.');
end
diffXY = norm(X-Y,'fro');
out.rel_diffXY = diffXY/max([1,norm(X,'fro'),norm(Y,'fro')]);
out.X = X;
out.value = sigma*norm(S*X-eye(n),'fro')^2/2 + norm(vec(X),1);

out.dualgap = dualgap;
end