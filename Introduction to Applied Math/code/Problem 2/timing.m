function u = timing(m,n,verbose)
%input m: size of temporal grid, n: size of spatial grid, verbose: name of
%the method.
%output u: results of the last layer.
%measuring time by tic toc.
    function B = cho(A)%cholesky decomposition
        le = length(A);
        for c=1:le 
            A(c:le, c) = A(c:le, c) - A(c:le, 1:c-1) * A(c, 1:c-1)';
            A(c, c) = sqrt(A(c, c));
            A(c+1:le, c) = A(c+1:le, c) / A(c, c); 
        end
        B = tril(A);
    end
k = 1/m;h = 1/n;
u = zeros((n-1)^2,1);
for i = 1:n-1
    for j = 1:n-1
        u((n-1)*(i-1)+j) = sin(pi*i/n)*sin(pi*j/n);
    end
end
%construct A.
for i = 1:(n-1)^2
    dig(i) = 1+4*k/h^2;
end
for i = 1:n^2-2*n
    subdig(i) = -k/h^2;
end
for i = 1:(n-1)*(n-2)
    subsubdig(i) = -k/h^2;
end
x_axis = [1:(n-1)^2,2:(n-1)^2,1:n^2-2*n,1:(n-1)*(n-2),n:(n-1)^2];
y_axis = [1:(n-1)^2,1:n^2-2*n,2:(n-1)^2,n:(n-1)^2,1:(n-1)*(n-2)];
value = [dig,subdig,subdig,subsubdig,subsubdig];
A = sparse(x_axis,y_axis,value);
for i = 1:n-2
    A((n-1)*i+1,(n-1)*i) = 0;
    A((n-1)*i,(n-1)*i+1) = 0;
end
if strcmp(verbose, 'squareroot')
    tic
    B = cho(A);
    for i = 1:m
        v = B\u;
        u = B'\v;
    end
    toc
end
if strcmp(verbose, 'gauss')
    tic
    L = tril(A);
    U = -triu(A,1);
    for i = 1:m
        u1 = u;
        init_err = norm(A*u-u);
        err = init_err;
        while err > 10^(-6)*init_err
            u1 = L\(U*u1 + u);
            err = norm(A*u1-u);
        end
        u = u1;
    end
    toc
end
if strcmp(verbose, 'gradient')
    tic
    for i = 1:m
        x0 = u;
        r0 = -A*u+u;
        init_err = norm(r0);
        err = init_err;
        p0 = r0;
        a0 = r0'*r0/(p0'*A*p0);
        x0 = x0 + a0*p0;
        r1 = r0 - a0*A*p0;
        while err > 10^(-6)*init_err
            b0 = (r1'*r1)/(r0'*r0);
            p0 = r1 + b0*p0;
            p1 = A*p0;
            a0 = (r1'*r1)/(p0'*p1);
            x0 = x0 + a0*p0;
            r0 = r1;
            r1 = r0 - a0*p1;
            err = norm(r1);
        end
        u = x0;
    end
    toc
end
%real solution.
X_real = zeros((n-1)^2,1);
for j = 1:(n-1)
    for p = 1:(n-1)
        X_real((n-1)*(j-1)+p) = sin(pi*p/n)*sin(pi*j/n)*exp(-2*pi^2);
    end
end
%err as output
err = norm(u-X_real);
end
