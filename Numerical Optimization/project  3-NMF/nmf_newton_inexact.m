function [x, infos] = nmf_newton_inexact(A, r, opts)
% Projected gradient descent for non-negative matrix factorization (NMF).
%
% The problem of interest is defined as
%
%           min 0.5 || A - BC ||_F^2,
%           where 
%           {A, B, C} > 0.
%
% Given a non-negative matrix A, factorized non-negative matrices {B, C} are calculated.
%
%
% Inputs:
%       A           : (m x n) non-negative matrix to factorize
%       r        : rank
%       in_options 
%           alg     : pgd: Projected gradient descent
%
%                   : direct_pgd: Projected gradient descent
%                       Reference:
%                           C.-J. Lin. 
%                           "Projected gradient methods for non-negative matrix factorization," 
%                           Neural Computation, vol. 19, pp.2756-2779, 2007.
%
% Output:
%       x           : non-negative matrix solution, i.e., x.B: (m x rank), x.C: (rank x n)
%       infos       : log information
%           epoch   : iteration nuber
%           cost    : objective function value
%           rel_cost: \|A-BC\|_F/\|A\|_F
%           rel_projnorm : projnorm/initial projnorm
%           time    : elapsed time
%           count   : number of inner iteration.
%
% Some codes are directly modifed from the original paper.

if nargin < 3
    opts = [];
end

if ~isfield(opts, 'max_epoch'); opts.max_epoch = 5e3; end
if ~isfield(opts, 'alpha'); opts.alpha = 0.1; end
if ~isfield(opts, 'tol_grad'); opts.tol_grad = 1e-4; end
if ~isfield(opts, 'tau'); opts.tau = 5; end

tic;

tau = opts.tau;
alpha = opts.alpha;
[m, n] = size(A);
epoch = 0;    
grad_count = 0;
B = rand(m, r); C = rand(r, n);

gradB = B*(C*C.') - A*C'; 
gradC = (B'*B)*C - B'*A;   
init_grad = norm([gradB; gradC'],'fro');


while 1
    projnorm = norm([gradB(gradB<0 | B>0); gradC(gradC<0 | C>0)],'fro');
    abs_cost = norm(A - B*C,'fro'); rel_cost = abs_cost/norm(A,'fro');
    fprintf('alg (Fast Newton): Epoch = %04d,\n relative cost = %.4e, relative projnorm = %.4e\n', epoch, rel_cost, projnorm/init_grad);
    if projnorm < opts.tol_grad*init_grad || epoch > opts.max_epoch
        break
    end
    
    C_old = C; BTB = B'*B;
    for i = 1:tau
        gradC_old = (B'*B)*C_old - B'*A; 
        I_plus = (abs(C_old) < eps) & (gradC_old > eps);
        U = (1-I_plus).*gradC_old; U = (1-I_plus).*(pinv(BTB)*U);
        C_old = max(C_old-alpha*U,0);
    end
    C = C_old; gradC = (B'*B)*C - B'*A; 
    B_old = B; CTC = C*C';
    for i = 1:tau
        gradB_old = B_old*(C*C') - A*C'; 
        I_plus = (abs(B_old) < eps) & (gradB_old > eps);
        U = (1-I_plus).*gradB_old; U = (1-I_plus).* (U*pinv(CTC));
        B_old = max(B_old-alpha*U,0);
    end
    B = B_old; gradB = B*(C*C') - A*C'; 
    epoch = epoch + 1;
    grad_count = grad_count + 1;
end
        

if epoch == opts.max_epoch
    fprintf('Max epoch reached.\n');
end
x.B = B;
x.C = C;
infos.epoch = epoch;
infos.time = toc;
infos.grad_count = grad_count;
infos.rel_cost = norm(A-B*C,'fro')/norm(A,'fro');
infos.rel_projnorm = projnorm/init_grad;
infos.projnorm = projnorm;

end

