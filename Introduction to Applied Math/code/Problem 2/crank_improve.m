function u = crank_improve(m,n)
%input m: size of temporal grid, n: size of spatial grid
%output u: results of the last layer.
%using conjugate gradient method.
%construct matrix.
k = 1/m;h = 1/n;
u = zeros((n-1)^2,1);
for i = 1:n-1
    for j = 1:n-1
        u((n-1)*(i-1)+j) = sin(pi*i/n)*sin(pi*j/n);
    end
end
for i = 1:(n-1)^2
    dig(i) = 1+2*k/h^2;
end
for i = 1:n^2-2*n
    subdig(i) = -k/2/h^2;
end
for i = 1:(n-1)*(n-2)
    subsubdig(i) = -k/2/h^2;
end
x_axis = [1:(n-1)^2,2:(n-1)^2,1:n^2-2*n,1:(n-1)*(n-2),n:(n-1)^2];
y_axis = [1:(n-1)^2,1:n^2-2*n,2:(n-1)^2,n:(n-1)^2,1:(n-1)*(n-2)];
value = [dig,subdig,subdig,subsubdig,subsubdig];
A = sparse(x_axis,y_axis,value);
for i = 1:n-2
    A((n-1)*i+1,(n-1)*i) = 0;
    A((n-1)*i,(n-1)*i+1) = 0;
end
for i = 1:(n-1)^2
    dig(i) = 1-2*k/h^2;
end
for i = 1:n^2-2*n
    subdig(i) = k/2/h^2;
end
for i = 1:(n-1)*(n-2)
    subsubdig(i) = k/2/h^2;
end
x_axis = [1:(n-1)^2,2:(n-1)^2,1:n^2-2*n,1:(n-1)*(n-2),n:(n-1)^2];
y_axis = [1:(n-1)^2,1:n^2-2*n,2:(n-1)^2,n:(n-1)^2,1:(n-1)*(n-2)];
value = [dig,subdig,subdig,subsubdig,subsubdig];
B = sparse(x_axis,y_axis,value);
for i = 1:n-2
    B((n-1)*i+1,(n-1)*i) = 0;
    B((n-1)*i,(n-1)*i+1) = 0;
end
for i = 1:m %adopt adjoint gradient descent.
    u = B*u;
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
end