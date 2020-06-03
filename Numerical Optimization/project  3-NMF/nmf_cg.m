function [x, infos] = nmf_cg(V, r, opts)
% Conjugate gradient descent for non-negative matrix factorization (NMF).
%
% The problem of interest is defined as
%
%           min || V - WH ||_F^2/2,
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
%           dfval   : difference of objective value
%           time    : elapsed time
%           grad_count : number of sampled data elements (gradient calculations)

%
% Some codes are directly modifed from the original paper.



if nargin < 3
    opts = [];
end
if ~isfield(opts, 'max_epoch'); opts.max_epoch = 1e4; end
if ~isfield(opts, 'tol_dfval'); opts.tol_dfval  = 1e-6; end
if ~isfield(opts, 'alpha'); opts.alpha = 1; end
if ~isfield(opts, 'tol_grad'); opts.tol_grad = 1e-4; end

Err = []; dErr = [];
tic;
[m, n] = size(V);
epoch = 0;    
grad_count = 0;
W = rand(m, r); H = rand(r, n);
dfval = inf; fval = inf;

gradW = W*(H*H') - V*H'; 
gradH = (W'*W)*H - W'*V;   
init_grad = norm([gradW; gradH'],'fro');

tolW = max(1e-3,opts.tol_grad)*init_grad; 
tolH = tolW;

stopping = 0;
while (dfval > opts.tol_dfval) && (epoch < opts.max_epoch) && ~stopping
    oldW = W; oldH = H;

    % stopping condition
    projnorm = norm([gradW(gradW<0 | W>0); gradH(gradH<0 | H>0)]);
    if projnorm < opts.tol_grad*init_grad
        stopping = 1;
    end
    
    [W, gradW, iterW] = nmf_subcg(V', H', W', tolW, 1000); 
    W = W'; 
    gradW = gradW';
    if iterW == 1
        %  the projected gradient method solves the subproblem stops without any iterations,
        tolW = 0.1 * tolW;
    end
    
    [H, gradH, iterH] = nmf_subcg(V, W, H, tolH, 1000);
    if iterH == 1
        tolH = 0.1 * tolH; 
    end
    W = max(W,0);H = max(H,0);
    
    fval0 = fval;
    fval = norm(V - W*H, 'fro')^2/2;
    Err = [Err fval];
    dfval = abs(fval0-fval)/max(fval,1);
    dErr = [Err fval];
    epoch = epoch + 1;
    grad_count = grad_count + m * n;
    
    fprintf('alg (%s): Epoch = %04d\n cost = %.8e, relative cost = %.4e, projnorm: %.4e\n', 'CG: LS', epoch, fval, dfval,projnorm);

end

if epoch == opts.max_epoch
    fprintf('Max epoch reached.\n');
end
x.W = W;
x.H = H;
infos.epoch = epoch;
infos.time = toc;
infos.grad_count = grad_count;
infos.dfval = dfval;
infos.fval = fval;

infos.rel_cost = norm(V - W*H,'fro')/norm(V,'fro');
infos.rel_projnorm = projnorm/init_grad;
end

function [H,grad,iter] = nmf_subcg(V,W,H,tol,maxiter)
[m,n] = size(H);
fun = @(H) sub_iter(V,W,H,m,n);
opts.maxIter = maxiter;
opts.tol = tol;
opts.debug = 0;
[H_0,out] = LS(fun, vec(H), opts);
H = reshape(H_0,m,n);
iter = out.iter;
grad = W'*W*H - W'*V;
end

function [f,g] = sub_iter(V,W,H_0,m,n)
H = reshape(H_0,m,n);
f = norm(V-W*H,'fro')^2/2;
g =  vec(W'*W*H - W'*V);
end