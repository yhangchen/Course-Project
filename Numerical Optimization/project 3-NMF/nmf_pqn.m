function [x, infos] = nmf_pqn(V, r, opts)
% quasi-Newton for non-negative matrix factorization (NMF).
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
if ~isfield(opts, 'alg'); opts.alg = 'pgd'; end
if ~isfield(opts, 'max_epoch'); opts.max_epoch = 5e3; end
if ~isfield(opts, 'alpha'); opts.alpha = 1; end
if ~isfield(opts, 'tol_grad'); opts.tol_grad = 1e-4; end

tic;
[m, n] = size(V);
epoch = 0;    
grad_count = 0;
W = rand(m, r); H = rand(r, n);

gradW = W*(H*H') - V*H'; 
gradH = (W'*W)*H - W'*V;   
init_grad = norm([gradW; gradH'],'fro');
tol_inner = 1e-8;

while epoch < opts.max_epoch
    % stopping condition
    projnorm = norm([gradW(gradW<0 | W>0); gradH(gradH<0 | H>0)]);
    if projnorm < opts.tol_grad*init_grad
        break
    end
    
    [W, gradW, iterW] = nmf_subpqn(V', H', W', tol_inner, 1000); 
    W = W'; 
    gradW = gradW';
%     iterW
%     if iterW == m
%         tolW
%         %  the projected gradient method solves the subproblem stops without any iterations,
%         tolW = 0.1 * tolW;
%     end

    [H, gradH, iterH] = nmf_subpqn(V, W, H, tol_inner, 1000);
%     if iterH == n
%         tolH = 0.1 * tolH; 
%     end
    gradW = W*(H*H') - V*H'; 
    grad_count = grad_count + iterW + iterH;

H = H + (H<eps) .* eps;
W = W + (W<eps) .* eps;
epoch = epoch + 1;
abs_cost = norm(V - W*H,'fro');
rel_cost = abs_cost/norm(V,'fro');
fprintf('alg (%s): Epoch = %04d,\n relative cost = %.4e, relative projnorm = %.4e\n', opts.alg, epoch, rel_cost, projnorm/init_grad);
end


if epoch == opts.max_epoch
    fprintf('Max epoch reached.\n');
end
x.W = W;
x.H = H;
infos.epoch = epoch;
infos.time = toc;
infos.grad_count = grad_count;
infos.rel_cost = rel_cost;
infos.rel_projnorm = projnorm/init_grad;
end

function [H, gradH, iterH] = nmf_subpqn(V,W,H,tol,maxIter)
    H = zeros(size(H)); gradH = zeros(size(H));
    iterH = 0;
    for i = 1:size(H,2)
        [h0, iterH0,  gradH0] = nnls_pqn(W, V(:,i), H(:,i), tol, maxIter);
        H(:,i) = h0;
        gradH(:,i) = gradH0;
        iterH = iterH + iterH0;
    end
end
