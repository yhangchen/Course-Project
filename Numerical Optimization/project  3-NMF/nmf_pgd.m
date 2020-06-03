function [x, infos] = nmf_pgd(V, r, opts)
% Projected gradient descent for non-negative matrix factorization (NMF).
%
% The problem of interest is defined as
%
%           min 0.5 || V - WH ||_F^2,
%           where 
%           {V, W, H} > 0.
%
% Given a non-negative matrix V, factorized non-negative matrices {W, H} are calculated.
%
%
% Inputs:
%       V           : (m x n) non-negative matrix to factorize
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
%       x           : non-negative matrix solution, i.e., x.W: (m x rank), x.H: (rank x n)
%       infos       : log information
%           epoch   : iteration nuber
%           cost    : objective function value
%           rel_cost: \|V-WH\|_F/\|V\|_F
%           rel_projnorm : projnorm/initial projnorm
%           time    : elapsed time
%           count   : number of inner iteration.
%
% Some codes are directly modifed from the original paper.



if nargin < 3
    opts = [];
end
if ~isfield(opts, 'alg'); opts.alg = 'pgd'; end
if ~isfield(opts, 'max_epoch'); opts.max_epoch = 5e3; end
if ~isfield(opts, 'tol_dfval'); opts.tol_dfval  = 1e-6; end
if ~isfield(opts, 'alpha'); opts.alpha = 1; end
if ~isfield(opts, 'tol_grad'); opts.tol_grad = 1e-4; end

if ~strcmp(opts.alg, 'pgd') && ~strcmp(opts.alg, 'direct_pgd')
    fprintf('Invalid algorithm: %s. Therfore, we use pgd.\n', opts.alg);
    opts.alg = 'pgd';
end

tic;
[m, n] = size(V);
epoch = 0;    
grad_count = 0;
W = rand(m, r); H = rand(r, n);

gradW = W*(H*H') - V*H'; 
gradH = (W'*W)*H - W'*V;   
init_grad = norm([gradW; gradH'],'fro');

if strcmp(opts.alg, 'pgd')
    tolW = max(1e-3,opts.tol_grad)*init_grad; 
    tolH = tolW;
else
    H = nmf_subpgd(V,W,H,0.001,1000);    
    obj = 0.5*(norm(V-W*H,'fro')^2);
    alpha = opts.alpha;
end

while epoch < opts.max_epoch
if strcmp(opts.alg, 'pgd')

    % stopping condition
    projnorm = norm([gradW(gradW<0 | W>0); gradH(gradH<0 | H>0)],'fro');
    if projnorm < opts.tol_grad*init_grad
        break
    end
    
    [W, gradW, iterW] = nmf_subpgd(V', H', W', tolW, 1000); 
    W = W'; 
    gradW = gradW';

    if iterW == 1
        %  the projected gradient method solves the subproblem stops without any iterations,
        tolW = 0.1 * tolW;
    end

    [H, gradH, iterH] = nmf_subpgd(V, W, H, tolH, 1000);
    if iterH == 1
        tolH = 0.1 * tolH; 
    end
    gradW = W*(H*H') - V*H'; 
    grad_count = grad_count + iterW + iterH;

elseif strcmp(opts.alg, 'direct_pgd')

    gradW = W*(H*H') - V*H';
    gradH = (W'*W)*H - W'*V;

    projnorm = norm([gradW(gradW<0 | W>0); gradH(gradH<0 | H>0)],'fro');  
    if projnorm < opts.tol_grad*init_grad
        fprintf('final grad norm %f\n', projnorm);
        break
    else
        Wn = max(W - alpha*gradW,0);    
        Hn = max(H - alpha*gradH,0);    
        newobj = 0.5*(norm(V-Wn*Hn,'fro')^2);
        inner_count = 0;
        if newobj-obj > 0.01*(sum(sum(gradW.*(Wn-W)))+ sum(sum(gradH.*(Hn-H))))
            % decrease stepsize    
            while 1
                alpha = alpha/10;
                Wn = max(W - alpha*gradW,0);    
                Hn = max(H - alpha*gradH,0);    
                newobj = 0.5*(norm(V-Wn*Hn,'fro')^2);
                inner_count = inner_count + 1;
                if newobj - obj <= 0.01*(sum(sum(gradW.*(Wn-W)))+ sum(sum(gradH.*(Hn-H))))
                    W = Wn; H = Hn;
                    obj = newobj;
                break;

                end
            end
        else 
            % increase stepsize
            while 1
                Wp = Wn; 
                Hp = Hn; 
                objp = newobj;
                alpha = alpha*10;
                Wn = max(W - alpha*gradW,0);    
                Hn = max(H - alpha*gradH,0);    
                newobj = 0.5*(norm(V-Wn*Hn,'fro')^2);
                inner_count = inner_count + 1;
                if (newobj - obj > 0.01*(sum(sum(gradW.*(Wn-W)))+ sum(sum(gradH.*(Hn-H))))) ...
                        || (isequal(Wn, Wp) && isequal(Hn, Hp))               
                    W = Wp; 
                    H = Hp;
                    obj = objp; 
                    alpha = alpha/10;
                    break;
                end
            end
        end
        grad_count = grad_count + inner_count;
    end

end
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
infos.rel_cost = norm(V - W*H,'fro')/norm(V,'fro');
infos.rel_projnorm = projnorm/init_grad;
infos.projnorm = projnorm;

end

function [H,grad,iter] = nmf_subpgd(V,W,H,tol,maxiter)
% implement the subprocess of PGD.
% H, grad: output solution and gradient
% iter: #iterations used
% V, W: constant matrices
% Hinit: initial solution
% tol: stopping tolerance
% maxiter: limit of iterations


% from the original code https://www.csie.ntu.edu.tw/~cjlin/papers/pgradnmf.pdf
WtV = W'*V;
WtW = W'*W; 

alpha = 1; beta = 0.1;
for iter=1:maxiter  
  grad = WtW*H - WtV;
  projgrad = norm(grad(grad<0|H>0));
  if projgrad < tol; break; end
  % search step size 
  for inner_iter=1:20
    Hn = max(H - alpha*grad, 0); d = Hn - H;
    gradd=sum(sum(grad.*d)); dQd = sum(sum((WtW*d).*d));
    suff_decr = 0.99*gradd + 0.5*dQd < 0;
    if inner_iter==1
      decr_alpha = ~suff_decr; Hp = H;
    end
    if decr_alpha
      if suff_decr
	H = Hn; break;
      else
	alpha = alpha * beta;
      end
    else
      if ~suff_decr | Hp == Hn
	H = Hp; break;
      else
	alpha = alpha/beta; Hp = Hn;
      end
    end
  end
end

if iter==maxiter
  fprintf('Max iter in nlssubprob\n');
end
end

