function [X, out] = spinv_alm(S, sigma, opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% not widely tested due to being inferior to spinv_admm.m

% (1) X^{k+1} := argmin_X f(X) +
%              g_sigma(Y^k)+<grad(g_sigma(Y^k)),X-Y^k>+||X-Y^k||_F^2/(2 mu)
%
% (2) Y^{k+1} := argmin_Y f(X^{k+1}) + <grad(f(X^{k+1})), Y-X^{k+1}>
%                         + ||Y -X^{k+1}||_F^2/(2 mu) + g_sigma(Y)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tuning mu
if ~isfield(opts,'mu_0'); opts.mu_0 = 0.1; end
if ~isfield(opts,'mu_N'); opts.mu_N = 40; end
if ~isfield(opts,'mu_ratio'); opts.mu_ratio = 3/4; end
if ~isfield(opts,'mu_lower'); opts.mu_lower = 1e-6; end
% stopping criterion
if ~isfield(opts,'tol_rel'); opts.tol_rel = eps; end 
if ~isfield(opts,'tol_gap'); opts.tol_gap = 1e-6; end
if ~isfield(opts,'maxIter'); opts.maxIter = 1e4; end
n = size(S,1);
% initial X: sparsest, initial Y: densest
X = zeros(n,n); Y = ones(n,n); Lambda = zeros(n,n);
mu = opts.mu_0;
dualgap = inf;
S2 = S*S; % store S^2
f1 = inf;
for iter = 1: opts.maxIter
    % update X
    Xp = X;
    B = mu*sigma/2*S2; A = B + eye(n);
    C = mu*(Lambda+sigma*S) + Y;
    X = sylvester(A,B,C);
    if sigma*norm(X,1) > sigma*norm(Y,1)-trace(Lambda'*(X-Y)) + norm(X-Y,'fro')^2/2/mu
        X = Y;
    end
    % update Y
    Yp = Y;
    Y0 = X - mu*Lambda;
    Y = sign(Y0).*max(abs(Y0)-mu,0);
    % update Lambda
    Lambda = Lambda - (X - Y)/mu/sigma;
    if mod(iter,opts.mu_N) == 0
        mu = max(mu * opts.mu_ratio, opts.mu_lower);
    end
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
out = [];
try
    chol(X+eps*eye(size(X)));
catch
    [V,D] = eig(X);
    d = diag(D); D = diag(max(d,0));
    X = (V*diag(d))*V';
end
out.X = X;
out.iter = iter; 
if iter == opts.maxIter
    fprintf('maximum number of iteration achieved.');
end
diffXY = norm(X-Y,'fro');
out.value = sigma*norm(S*X-eye(n),'fro')^2/2 + norm(vec(X),1);
out.rel_diffXY = diffXY/max([1,norm(X,'fro'),norm(Y,'fro')]);
out.gap = dualgap;
end
