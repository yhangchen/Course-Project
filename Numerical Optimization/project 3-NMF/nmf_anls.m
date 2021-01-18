function [x, infos] = nmf_anls(V, r, opts)
% Alternative non-negative least squares (ANLS) for non-negative matrix factorization (NMF).
%
% The problem of interest is defined as
%
%           min 1/2 || V - WH ||_F^2,
%           where 
%           {V, W, H} > 0.
%
% Given a non-negative matrix V, factorized non-negative matrices {W, H} are calculated.
%
%
% Inputs:
%       V           : (m x n) non-negative matrix to factorize
%       rank        : rank
%       in_opts  : opts
%
%
% References:
%       Jingu Kim, Yunlong He, and Haesun Park,
%       "Algorithms for Nonnegative Matrix and Tensor Factorizations: A Unified View 
%       Based on Block Coordinate Descent Framework,"
%       Journal of Global Optimization, 58(2), pp. 285-319, 2014.
%
%       Jingu Kim and Haesun Park.
%       "Fast Nonnegative Matrix Factorization: An Active-set-like Method and Comparisons,"
%       SIAM Journal on Scientific Computing (SISC), 33(6), pp. 3261-3281, 2011.
%
%
% Output:
%       x           : non-negative matrix solution, i.e., x.W: (m x rank), x.H: (rank x n)
%       infos       : log information
%           epoch   : iteration nuber
%           cost    : objective function value
%           optgap  : optimality gap
%           time    : elapsed time
%           grad_calc_count : number of sampled data elements (gradient calculations)

    [m, n] = size(V);
 
    % set local opts
    if nargin < 3
        opts = [];
    end
    % merge opts
    if ~isfield(opts, 'alg'); opts.alg = 'anls_asgivens'; end
    if ~isfield(opts, 'max_epoch'); opts.max_epoch = 5e3; end
    if ~isfield(opts, 'alpha'); opts.alpha = 1; end
    if ~isfield(opts, 'tol_grad'); opts.tol_grad = 1e-4; end

% set paramters
    if ~strcmp(opts.alg, 'anls_asgroup') && ~strcmp(opts.alg, 'anls_asgivens') ...
            && ~strcmp(opts.alg, 'anls_bpp') 
        fprintf('Invalid algorithm: %s. Therfore, we use anls_asgroup (i.e., ANLS with Active Set Method and Column Grouping).\n', opts.alg);
        opts.alg = 'anls_asgroup';
    end
    
    W = rand(m, r); H = rand(r, n);

    gradW = W*(H*H') - V*H'; 
    gradH = (W'*W)*H - W'*V;   
    init_grad = norm([gradW; gradH'],'fro');

            

    tic;
    
    epoch = 0;
    
    % main loop
    while (epoch < opts.max_epoch)
        gradW = W*(H*H') - V*H'; gradH = (W'*W)*H - W'*V;
        projnorm = norm([gradW(gradW<0 | W>0); gradH(gradH<0 | H>0)],'fro');
        if projnorm < opts.tol_grad*init_grad
            break
        end
        if strcmp(opts.alg, 'anls_asgroup')
            ow = 0;
            H = nnlsm_activeset(W'*W, W'*V, ow, 1, H);
            W = nnlsm_activeset(H*H', H*V',ow, 1, W');
            W = W';
            
        elseif strcmp(opts.alg, 'anls_asgivens')
            ow = 0;
            WtV = W' * V;
            for i=1:size(H,2)
                H(:,i) = nnls1_asgivens(W'*W, WtV(:,i), ow, 1, H(:,i));
            end

            HAt = H*V';
            Wt = W';
            for i=1:size(W,1)
                Wt(:,i) = nnls1_asgivens(H*H', HAt(:,i), ow, 1, Wt(:,i));
            end
            W = Wt';
            
        elseif strcmp(opts.alg, 'anls_bpp')
            H = nnlsm_blockpivot(W'*W, W'*V, 1, H);
            W = nnlsm_blockpivot(H*H', H*V', 1, W');
            W = W';

        end
        

        epoch = epoch + 1;         
        
        abs_cost = norm(V - W*H,'fro');
        rel_cost = abs_cost/norm(V,'fro');
        fprintf('alg (%s): Epoch = %04d,\n relative cost = %.4e, relative projnorm = %.4e\n', opts.alg, epoch, rel_cost, projnorm/init_grad);

    end
    
    x.W = W;
    x.H = H;
infos.epoch = epoch; infos.time = toc; infos.projnorm = projnorm; 
infos.rel_cost = norm(V - W*H,'fro')/norm(V,'fro');
infos.rel_projnorm = projnorm/init_grad;
infos.projnorm = projnorm;
end
