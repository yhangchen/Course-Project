function [x, infos] = nmf_mu(V, r, opts)
% Multiplicative upates (MU) for non-negative matrix factorization (NMF).
%
% The problem of interest is defined as
%
%           min f(V, W, H),
%           where 
%           {V, W, H} >= 0.
%
% Given a non-negative matrix V, factorized non-negative matrices {W, H} are calculated.
%
%
% Inputs:
%       V           : (m x n) non-negative matrix to factorize
%       r        : rank
%       opts 
%             metric: 
%                 'fro' f(V,W,H) = ||V-W*H||_F^2 / 2
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
if ~isfield(opts, 'alg'); opts.alg = 'mu'; end
if ~isfield(opts, 'metric'); opts.metric = 'fro'; end
if ~isfield(opts, 'max_epoch'); opts.max_epoch = 5e3; end
if ~isfield(opts, 'tol_grad'); opts.tol_grad = 1e-4; end


[m, n] = size(V);
epoch = 0;    
count = 0;
W = rand(m, r); H = rand(r, n);


gradW = W*(H*H') - V*H'; 
gradH = (W'*W)*H - W'*V;   
init_grad = norm([gradW; gradH'],'fro');

tic
while (epoch < opts.max_epoch)
% if strcmp(opts.metric, 'fro')
    % update H
    H = H .* (W' * V) ./ (W' * W * H);
    H = H + (H<eps) .* eps;
    % update W
    W = W .* (V * H') ./ (W * (H * H'));
    W = W + (W<eps) .* eps;
    % calulate error
    gradW = W*(H*H') - V*H';
    gradH = (W'*W)*H - W'*V;
    projnorm = norm([gradW(gradW<0 | W>0); gradH(gradH<0 | H>0)],'fro');

% elseif strcmp(opts.metric, 'KL')
%     % update W
%     W = W .* ((V./(W*H + eps))*H')./(ones(m,1)*sum(H'));
%     % update H
%     H = H .* (W'*(V./(W*H + eps)))./(sum(W)'*ones(1,n));
%     % calulate error
%     V1 = W * H;
%     V1 = V1 + (V1<eps) .* eps;
%     temp = V.*log(V./V1);
%     temp(temp ~= temp) = 0;
%     fval0 = fval;
%     fval = sum(sum(temp - V + V1));
% else
%     error('undefined metric')
% end
abs_cost = norm(V - W*H,'fro');
rel_cost = abs_cost/norm(V,'fro');
epoch = epoch + 1;
count = count + 1;
    if projnorm < opts.tol_grad*init_grad
        break
    end
fprintf('metric (%s): Epoch = %04d,\n relative cost = %.4e, relative projnorm: %.4e\n', opts.metric, epoch, rel_cost, projnorm/init_grad);
end
if epoch == opts.max_epoch
    fprintf('Max epoch reached.\n');
end



x.W = W;
x.H = H;
infos.epoch = epoch;
infos.time = toc;
infos.count = count;
infos.rel_cost = norm(V - W*H,'fro')/norm(V,'fro');
infos.rel_projnorm = projnorm/init_grad;
infos.projnorm = projnorm;

end

