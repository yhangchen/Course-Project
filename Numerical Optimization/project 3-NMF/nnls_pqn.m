function [x, epoch, grad] = nnls_pqn(A, b, x0, tol, maxIter)
%
% function out = pqn(A, b, x0, options)
% Solve a nonnegative least squares problem, min ||Ax - b||_2^2 /2, s.t. x \ge 0
% Projected Quasi-Newton Based.
%
    if nargin < 5
        maxIter = 1000;
    end
    if nargin < 4
        tol = 1e-3;
    end
    
    tic
    x = x0; grad = A'*(A*x-b); obj = norm(A*x-b,'fro')^2/2;
    epoch = 0; inv_H = eye(length(x)); init_grad = norm(grad(grad<0|x>0));
    
    % Main loop
    while 1
        epoch = epoch + 1;
        if epoch ==  maxIter; break; end
        x_old = x; 
        I_p = x == 0 & grad > 0;
        grad(I_p) = 0;
        direc = -inv_H * grad;
        ns = x == 0 & direc < 0;
        direc(I_p) = 0;               % d
        direc(ns) = 0;               % \bar{d}

        Ad = A * direc;
        step = -Ad' * ((A * x) - b) / (Ad' * Ad);
        
        x = max(x + step * direc,0);
        dx = x - x_old;
        dx(I_p) = 0;
        
        grad = A'*(A*x-b); obj = norm(A*x-b,'fro')^2/2;
        grad_norm = norm(grad(grad<0|x>0));
        
        if grad_norm < tol * init_grad; break; end
        
        % BFGS
        if epoch > 4 
            Au = A * dx;
            norm_Au = norm(Au)^2;
            AtAu = A' * Au;
            DAtAu = inv_H * AtAu;
            inv_H = inv_H + ((1 + AtAu' * DAtAu / norm_Au) * (dx * dx') - (DAtAu * dx' + dx * DAtAu')) / norm_Au;
        end
    end
end
 