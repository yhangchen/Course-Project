function [u,err] = implicit(m, n)
%input m: size of temporal grid, n: size of spatial grid
%output u: results of the last layer, err: L2 error to the real solution.
%construct A.
k = 1/m;h = 1/n;
u = zeros((n-1)^2,1);
for i = 1:n-1
    for j = 1:n-1
        u((n-1)*(i-1)+j) = sin(pi*i/n)*sin(pi*j/n);
    end
end
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
for i = 1:m
    u = A\u;
end
X_real = zeros((n-1)^2,1);
for j = 1:(n-1)
    for p = 1:(n-1)
        X_real((n-1)*(j-1)+p) = sin(pi*p/n)*sin(pi*j/n)*exp(-2*pi^2);
    end
end
err = norm(u-X_real);
end
