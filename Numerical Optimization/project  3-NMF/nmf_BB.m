function [x, infos] = nmf_BB(V, r, opts)
% Projected gradient descent with BB stepsize for non-negative matrix factorization (NMF).
%
% The problem of interest is defined as
%
%           min || V - WH ||_F^2,
%           where 
%           {V, W, H} > 0.
%
% Given a non-negative matrix V, factorized non-negative matrices {W, H} are calculated.
%
%
% Inputs:
%       V           : (m x n) non-negative matrix to factorize
%       r        : rank
%
% Output:
%       x           : non-negative matrix solution, i.e., x.W: (m x rank), x.H: (rank x n)
%       infos       : log information
%           epoch   : iteration nuber
%           cost    : objective function value
%           rel_cost: \|V-WH\|_F/\|V\|_F
%           rel_projnorm : projnorm/initial projnorm
%           time    : elapsed time
%           count   : number of inner iteration.

if nargin < 3
    opts = [];
end
if ~isfield(opts,'max_epoch'); opts.max_epoch = 5e3; end
if ~isfield(opts, 'tol_grad'); opts.tol_grad = 1e-4; end

W = rand(size(V,1),r); H = rand(r,size(V,2));

tic;

epoch=0; count = 0;

gradW = W*(H*H') - V*H'; 
gradH = (W'*W)*H - W'*V;

init_grad = norm([gradW; gradH'],'fro');
tolW = max(1e-3, opts.tol_grad)*init_grad; tolH = tolW;

while (epoch < opts.max_epoch)

    projnorm = norm([gradW(gradW<0 | W>0); gradH(gradH<0 | H>0)]);
    if projnorm < opts.tol_grad*init_grad
        break
    end
    
    [W, gradW, iterW] = nmf_subBB(V',H',W',tolW,1000); 
    W = W'; gradW = gradW';
    if iterW == 1; tolW = 0.1 * tolW; end 

    [H, gradH, iterH] = nmf_subBB(V,W,H,tolH,1000);
    if iterH == 1; tolH = 0.1 * tolH; end
    
    epoch = epoch + 1;
    count = count + iterW + iterH;
    
    gradW = W*(H*H') - V*H'; 
    abs_cost = norm(V - W*H,'fro');
    rel_cost = abs_cost/norm(V,'fro');
    fprintf('alg (BB): Epoch = %04d\n, relative cost = %.4e, relative projnorm = %.4e\n', epoch, rel_cost, projnorm/init_grad);

end

x.W = W; x.H = H;
infos.epoch = epoch; infos.time = toc; infos.projnorm = projnorm; 
infos.count = count;
infos.rel_cost = rel_cost;
infos.rel_projnorm = projnorm/init_grad;
infos.projnorm = projnorm;

end

function [x, gradx, iter] = nmf_subBB(b, A, x0, tol, maxiter)
% Projected Barzilai-Borwein method for the 
% Nonnegative Least Squares Problem: min 0.5 \|Ax-b\|_{F}^2 subject to x>=0
% using a nonmonotone line search.
%
%%%%%%%%%%%%%%%%%%%%%%%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%
% B: m-by-n matrix; A: m-by-r matrix; 
% x0: random r-by-n initial solution; 
% tol: tolerance for stopping; 
% maxiter: maximum number of iterations.
%%%%%%%%%%%%%%%%%%%%%%%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%
% x: output; 
% gradx: gradient of objective function;
% iter: number of iterations that has been used.

M=10;  %mm nonmonotone line search parameter

gamma=1e-4;

x = x0; AtA=A'*A; Atb=A'*b;
gradx = AtA*x - Atb; 

f0 = (sum(sum(b.*b)) - 2*sum(sum(x.*Atb)) + sum(sum((AtA*x).*x))) / 2;

fs = [];
for iter=1:maxiter
    projnorm = norm(gradx(gradx<0 | x>0));

    if projnorm < tol; break; end

    if iter==1; fs(iter)=f0;
    else; fs(iter)=fn;
    end

    maxima = max(fs(max(0,end-M)+1:end));

    if iter==1; lambda=1/max(max(abs(gradx))); end
    x1 = max(x - lambda*gradx, 0); 
    dx = x1 - x;
    decres = sum(sum(dx.*gradx));
    dAd = sum(sum((AtA*dx).*dx));
    alpha=1;
    fn = fs(iter) + alpha*decres + 0.5*alpha^2*dAd; 
    while (fn > maxima + alpha*gamma*decres)
    % backtracking Line Search
        alpha = alpha/4;
        fn = fs(iter) + alpha*decres + 0.5*alpha^2*dAd;
    end
    x1 = x+alpha*dx;
    % Compute the BB steplength   
    gradx1 = AtA*x1-Atb;
    s1 = x1-x;  y1 = gradx1-gradx;
    sTs=sum(sum(s1.*s1));
    sTy=sum(sum(s1.*y1));
    if sTy <= 0
      lambda=1/eps;
    else
      lambda=min(1/eps,max(eps,sTs/sTy));
    end
    x = x1; 
    gradx = gradx1;

end
end