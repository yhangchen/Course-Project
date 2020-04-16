function [x,out] = l1l1_lbreg_bb(x0, A, b, mu, opts)
    % implement Linearized Bregman with BB from the paper
    % A linearized Bregman algorithm for decentralized basis pursuit
    
    [m0, n0] = size(A);
    A = [A mu*eye(m0)];
    b = mu*b;

    if isfield(opts,'tol'),    tol = opts.tol;     else; tol = 1e-12;   end
    if isfield(opts,'maxit'),  maxit = opts.maxit; else; maxit = 1e4;  end
    if isfield(opts,'x_real'),  x_real = opts.x_real; else; x_real = []; end
    if isfield(opts,'alpha'),  alpha = opts.alpha; else; alpha = 5*norm(x_real,inf); end
    if isfield(opts,'stepsize'),  stepsize0 = opts.stepsize;  else; stepsize0 = 2/alpha/normest(A*A',1e-1); end

    [m,~] = size(A);
    y = zeros(m,1);
    res = b;
    tolb = tol*norm(b);
    soft = @(z) sign(z).*max(abs(z)-1,0);
    for iter = 1:maxit
        if iter > 1
            dy = y - yp;
            dres = resp-res;
            stepsize = dy'*dy / (dy'*dres);
            if isinf(stepsize)
                stepsize = stepsize0; 
            end 
        else
            stepsize = stepsize0 + 1/max(abs(b'*A));
        end
        yp = y;
        y = y + stepsize * res; 
        x = alpha * soft(A'*y);
        resp = res; 
        res = b - A*x; 
        if iter > 1 && norm(res) < tolb
            break; 
        end
    end
    out.iter = iter;
    x = x(1:n0)/mu;
    if iter == maxit
        fprintf('have achieved maximum number of iteraton\n')
    end
end