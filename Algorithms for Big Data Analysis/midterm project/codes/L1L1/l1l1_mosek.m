function [x,out]=l1l1_mosek(x0, A, b, mu, opts)
    tic
	[m, n] = size(A);
    f1 = ones(2*n,1); f2 = ones(2*m,1);
    f = [mu*f1;f2];
    A1 = [A -A -eye(m) eye(m)];
    l = zeros(2*n+2*m,1);
    x0 = linprog(f,[],[],A1,b,l); % linprog in mosek
    
    % another formuation is slower.
    % f1 = ones(2*n,1); f2 = ones(m,1);
    % f = [mu*f1;f2];
    % A0 = [A -A -eye(m)];
    % A00 = [-A A -eye(m)];
    % A1 = [A0;A00];
    % b1 = [b;-b];
    % l = zeros(2*n+m,1);
    % x0 = linprog(f,A1,b1,[],[],l);

    time = toc;
	x = x0(1:n)-x0(n+1:2*n);
	out = [];
	out.time = time;
	out.val = norm(A*x-b,1)+mu*norm(x,1);
end


